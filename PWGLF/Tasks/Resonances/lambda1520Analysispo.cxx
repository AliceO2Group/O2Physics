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

/// \file Lambda1520Analysispo.cxx
/// \brief Task for Lambda(1520) resonance reconstruction via proton-kaon invariant mass analysis
///
/// \author Yash Patley <yash.patley@cern.ch>
/// \author Durgesh Bhatt <durgesh.bhatt@cern.ch>

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/Core/RecoDecay.h"

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

#include <cmath>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct Lambda1520Analysispo {

  SliceCache sliceCache;

  // ── Named PDG codes (fixes pdg/explicit-code and magic-number linter errors) ─
  static constexpr int kPdgProton{2212};
  static constexpr int kPdgKaon{321};
  static constexpr int kPdgLambda0{3122};
  static constexpr int kPdgXiMinus{3312};
  static constexpr int kPdgXi0{3322};
  static constexpr int kPdgOmegaMinus{3334};

  // Preslice helpers: allow fast lookup of tracks belonging to a collision
  Preslice<aod::ResoTracks> tracksPerResonanceCollision = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> tracksPerStandardCollision = aod::track::collisionId;

  // Pointer to MC parent particle table (used only in MC processing)
  aod::ResoMCParents const* mcResonanceParentTable = nullptr;

  // ── Event-level configurables ────────────────────────────────────────────
  Configurable<bool> applyOccupancyInTimeRangeCut{"applyOccupancyInTimeRangeCut", false, "If true, apply a cut on the number of tracks in a time window around the collision (occupancy cut)"};

  // ── Histogram binning configurables ─────────────────────────────────────
  Configurable<int> numberOfPtBins{"numberOfPtBins", 100, "Number of bins along the transverse momentum (pT) axis"};
  Configurable<int> numberOfInvMassBins{"numberOfInvMassBins", 120, "Number of bins along the invariant mass axis"};

  // ── Physics configurables ────────────────────────────────────────────────
  Configurable<int> pdgCodeLambda1520{"pdgCodeLambda1520", 3124, "PDG code of the Lambda(1520) resonance "};
  Configurable<bool> enableRotationalBackground{"enableRotationalBackground", true, "If true, compute rotational background (kaon phi rotated by ~pi) to estimate combinatorial background"};

  // ── Track quality cuts ───────────────────────────────────────────────────
  Configurable<float> minTrackPt{"minTrackPt", 0.15f, "Minimum transverse momentum pT of a track [GeV/c]"};
  Configurable<float> minTrackMomentum{"minTrackMomentum", 0.f, "Minimum total momentum p of a track [GeV/c]"};
  Configurable<float> minPseudorapidity{"minPseudorapidity", -0.8f, "Minimum pseudorapidity eta (detector acceptance)"};
  Configurable<float> maxPseudorapidity{"maxPseudorapidity", 0.8f, "Maximum pseudorapidity eta (detector acceptance)"};
  Configurable<float> maxDCAz{"maxDCAz", 1.0f, "Maximum allowed distance of closest approach along z-axis (DCAz) [cm]."};
  Configurable<float> minPairRapidity{"minPairRapidity", -0.5f, "Minimum rapidity y of the reconstructed proton-kaon pair"};
  Configurable<float> maxPairRapidity{"maxPairRapidity", 0.5f, "Maximum rapidity y of the reconstructed proton-kaon pair"};

  // ── TPC cluster quality cuts ─────────────────────────────────────────────
  Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 70, "Minimum number of TPC crossed pad rows (track quality)"};
  Configurable<bool> applyCrossedRowsCut{"applyCrossedRowsCut", false, "If true, require at least minTPCCrossedRows crossed rows in the TPC"};
  Configurable<int> minTPCClustersFound{"minTPCClustersFound", 70, "Minimum number of TPC clusters found on the track"};
  Configurable<bool> applyTPCClustersCut{"applyTPCClustersCut", false, "If true, require at least minTPCClustersFound TPC clusters"};

  // ── pT-dependent DCAxy cuts for Protons ─────────────────────────────────
  Configurable<std::vector<float>> protonDCAPtBinEdges{"protonDCAPtBinEdges", {0.0f, 0.5f, 1.0f, 2.0f, 3.0f, 5.0f, 1000.0f}, "Proton pT bin edges [GeV/c] for the pT-dependent DCAxy selection"};
  Configurable<std::vector<float>> protonMaxDCAxyPerPtBin{"protonMaxDCAxyPerPtBin", {0.020f, 0.015f, 0.010f, 0.007f, 0.005f, 0.004f}, "Maximum |DCAxy| [cm] for protons in each pT bin defined by protonDCAPtBinEdges"};

  // ── pT-dependent DCAxy cuts for Kaons ───────────────────────────────────
  Configurable<std::vector<float>> kaonDCAPtBinEdges{"kaonDCAPtBinEdges", {0.0f, 0.3f, 0.6f, 1.0f, 2.0f, 1000.0f}, "Kaon pT bin edges [GeV/c] for the pT-dependent DCAxy selection"};
  Configurable<std::vector<float>> kaonMaxDCAxyPerPtBin{"kaonMaxDCAxyPerPtBin", {0.025f, 0.018f, 0.012f, 0.008f, 0.004f}, "Maximum |DCAxy| [cm] for kaons in each pT bin defined by kaonDCAPtBinEdges"};

  // ── Analysis mode switches ───────────────────────────────────────────────
  Configurable<bool> runQualityChecksOnly{"runQualityChecksOnly", false, "If true, only fill QA histograms and skip invariant mass computation"};
  Configurable<bool> applyDeepAngleCut{"applyDeepAngleCut", false, "If true, reject proton-kaon pairs with very small opening angle (removes split-track background)"};
  Configurable<double> deepAngleCutValue{"deepAngleCutValue", 0.04, "Minimum allowed opening angle [rad] between proton and kaon (used if applyDeepAngleCut = true)"};
  Configurable<bool> applyKinematicPairCuts{"applyKinematicPairCuts", false, "If true, apply additional kinematic cuts on the p-K opening angle"};

  // ── Global track selection flags ─────────────────────────────────────────
  Configurable<bool> requirePrimaryTrack{"requirePrimaryTrack", true,
                                         "Require track to pass the 'isPrimaryTrack' flag (kGoldenChi2 | kDCAxy | kDCAz)"};
  Configurable<bool> requireGlobalTrackNoDCA{"requireGlobalTrackNoDCA", true,
                                             "Require track to pass 'isGlobalTrackWoDCA' (quality cuts without DCA requirement)"};
  Configurable<bool> requirePVContributor{"requirePVContributor", true,
                                          "Require track to be a contributor to the primary vertex reconstruction"};

  // ── PID configurables ────────────────────────────────────────────────────
  Configurable<bool> requireTOFForProton{"requireTOFForProton", false,
                                         "If true, only accept proton candidates that have a TOF measurement"};
  Configurable<bool> requireTOFForKaon{"requireTOFForKaon", false,
                                       "If true, only accept kaon candidates that have a TOF measurement"};
  Configurable<bool> useTPCOnlyPID{"useTPCOnlyPID", false,
                                   "If true, use only TPC nSigma for PID (ignore TOF even if available)"};

  Configurable<float> tpcNSigmaVetoOtherSpecies{"tpcNSigmaVetoOtherSpecies", 3.0f,
                                                "Reject track if its TPC nSigma for a different species is below this threshold (avoids misID)"};
  Configurable<float> tpcNSigmaVetoPion{"tpcNSigmaVetoPion", 3.0f,
                                        "TPC nSigma threshold below which a track is vetoed as a pion"};
  Configurable<float> tpcNSigmaVetoKaonForProton{"tpcNSigmaVetoKaonForProton", 3.0f,
                                                 "TPC nSigma threshold: veto proton candidates that look like kaons"};
  Configurable<float> tpcNSigmaVetoPionForKaon{"tpcNSigmaVetoPionForKaon", 3.0f,
                                               "TPC nSigma threshold: veto kaon candidates that look like pions"};
  Configurable<float> tpcNSigmaVetoProtonForKaon{"tpcNSigmaVetoProtonForKaon", 3.0f,
                                                 "TPC nSigma threshold: veto kaon candidates that look like protons"};

  Configurable<float> minTPCNSigmaKaon{"minTPCNSigmaKaon", -6.0f, "Minimum (most negative) TPC nSigma for kaon"};
  Configurable<float> minTPCNSigmaProton{"minTPCNSigmaProton", -6.0f, "Minimum (most negative) TPC nSigma for proton"};
  Configurable<float> minTOFNSigmaKaon{"minTOFNSigmaKaon", -6.0f, "Minimum (most negative) TOF nSigma for kaon"};
  Configurable<float> minTOFNSigmaProton{"minTOFNSigmaProton", -6.0f, "Minimum (most negative) TOF nSigma for proton"};
  Configurable<float> minCombinedTPCTOFNSigmaKaon{"minCombinedTPCTOFNSigmaKaon", -6.0f,
                                                  "Minimum combined TPC+TOF nSigma for kaon (used in combined PID mode)"};
  Configurable<float> minCombinedTPCTOFNSigmaProton{"minCombinedTPCTOFNSigmaProton", -6.0f,
                                                    "Minimum combined TPC+TOF nSigma for proton (used in combined PID mode)"};

  Configurable<float> tpcNSigmaVetoThreshold{"tpcNSigmaVetoThreshold", 3.0f,
                                             "General TPC nSigma veto cut to reject misidentified particles when TPC bands overlap"};
  Configurable<float> tofNSigmaVetoThreshold{"tofNSigmaVetoThreshold", 3.0f,
                                             "General TOF nSigma veto cut to reject misidentified particles"};

  // ── Proton PID momentum-dependent TPC cuts ───────────────────────────────
  Configurable<double> maxTPCNSigmaProton{"maxTPCNSigmaProton", 3.0,
                                          "Maximum |TPC nSigma| for proton identification (symmetric cut)"};
  Configurable<double> combinedNSigmaCutProton{"combinedNSigmaCutProton", 3.0,
                                               "Cut on sqrt(nSigmaTPC^2 + nSigmaTOF^2) for proton. Negative value switches to asymmetric mode."};
  Configurable<std::vector<float>> protonTPCPIDMomentumBins{"protonTPCPIDMomentumBins",
                                                            {0, 0.5, 0.7, 0.8},
                                                            "Momentum p bin edges [GeV/c] for momentum-dependent TPC PID cuts on protons"};
  Configurable<std::vector<float>> protonTPCNSigmaCutPerBin{"protonTPCNSigmaCutPerBin",
                                                            {5., 3.5, 2.5},
                                                            "Maximum TPC nSigma for proton in each momentum bin (tighter at higher p)"};
  Configurable<std::vector<float>> protonTOFPIDMomentumBins{"protonTOFPIDMomentumBins",
                                                            {0., 999.},
                                                            "Momentum p bin edges [GeV/c] for momentum-dependent TOF PID cuts on protons"};
  Configurable<std::vector<float>> protonTOFNSigmaCutPerBin{"protonTOFNSigmaCutPerBin",
                                                            {3.0},
                                                            "Maximum TOF nSigma for proton in each momentum bin"};

  // ── Kaon PID momentum-dependent TPC cuts ────────────────────────────────
  Configurable<double> maxTPCNSigmaKaon{"maxTPCNSigmaKaon", 3.0,
                                        "Maximum |TPC nSigma| for kaon identification (symmetric cut)"};
  Configurable<double> combinedNSigmaCutKaon{"combinedNSigmaCutKaon", 3.0,
                                             "Cut on sqrt(nSigmaTPC^2 + nSigmaTOF^2) for kaon. Negative value switches to asymmetric mode."};
  Configurable<std::vector<float>> kaonTPCPIDMomentumBins{"kaonTPCPIDMomentumBins",
                                                          {0., 0.25, 0.3, 0.45},
                                                          "Momentum p bin edges [GeV/c] for momentum-dependent TPC PID cuts on kaons"};
  Configurable<std::vector<float>> kaonTPCNSigmaCutPerBin{"kaonTPCNSigmaCutPerBin",
                                                          {6, 3.5, 2.5},
                                                          "Maximum TPC nSigma for kaon in each momentum bin"};
  Configurable<std::vector<float>> kaonTOFPIDMomentumBins{"kaonTOFPIDMomentumBins",
                                                          {0., 999.},
                                                          "Momentum p bin edges [GeV/c] for momentum-dependent TOF PID cuts on kaons"};
  Configurable<std::vector<float>> kaonTOFNSigmaCutPerBin{"kaonTOFNSigmaCutPerBin",
                                                          {3.0},
                                                          "Maximum TOF nSigma for kaon in each momentum bin"};

  // ── Event mixing configurables ───────────────────────────────────────────
  Configurable<int> numberOfEventsToMix{"numberOfEventsToMix", 20,
                                        "Number of events to mix with each signal event for background estimation"};
  ConfigurableAxis dcaZMixingBins{"dcaZMixingBins",
                                  {VARIABLE_WIDTH, -1.2f, -1.0f, -0.9f, -0.8f, -0.7f, -0.6f, -0.5f, -0.4f, -0.3f,
                                   -0.2f, -0.1f, 0.f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.2f},
                                  "DCAz vertex bins used for event mixing pairing"};
  ConfigurableAxis vertexZMixingBins{"vertexZMixingBins",
                                     {VARIABLE_WIDTH, -10.f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f,
                                      0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f},
                                     "Primary vertex z-position bins used for event mixing pairing"};
  ConfigurableAxis centralityMixingBins{"centralityMixingBins",
                                        {VARIABLE_WIDTH, 0.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f, 100.f, 200.f},
                                        "Centrality (FT0 %) bins used for event mixing pairing"};
  ConfigurableAxis eventPlaneMixingBins{"eventPlaneMixingBins",
                                        {VARIABLE_WIDTH, -1.5708f, -1.25664f, -0.942478f, -0.628319f, 0.f,
                                         0.628319f, 0.942478f, 1.25664f, 1.5708f},
                                        "Event plane angle bins used for event mixing (EP-dependent analysis)"};
  ConfigurableAxis occupancyBins{"occupancyBins",
                                 {VARIABLE_WIDTH, 0.0, 100, 500, 600, 1000, 1100, 1500, 1600, 2000, 2100, 2500,
                                  2600, 3000, 3100, 3500, 3600, 4000, 4100, 4500, 4600, 5000, 5100, 9999},
                                 "Track occupancy bins in the time range around the collision"};

  // ── Rotational background configurables ─────────────────────────────────
  Configurable<int> numberOfRotations{"numberOfRotations", 10,
                                      "How many times to rotate the kaon phi for the rotational background estimate"};
  Configurable<float> rotationAngleWindow{"rotationAngleWindow", 6.f,
                                          "The kaon is rotated by angles near PI, within a window of PI/rotationAngleWindow"};

  // ── MC event selection flags ─────────────────────────────────────────────
  Configurable<bool> mcRequireAfterAllCuts{"mcRequireAfterAllCuts", false,
                                           "MC event selection: require isInAfterAllCuts flag"};
  Configurable<bool> mcRequireINELgt0{"mcRequireINELgt0", false,
                                      "MC event selection: require at least 1 charged particle in |eta|<1 (INEL>0)"};
  Configurable<bool> mcRequireSel8{"mcRequireSel8", false,
                                   "MC event selection: require the standard Sel8 event selection"};
  Configurable<bool> mcRequireVtxWithin10cm{"mcRequireVtxWithin10cm", false,
                                            "MC event selection: require primary vertex |z| < 10 cm"};
  Configurable<bool> mcRequireTriggerTVX{"mcRequireTriggerTVX", false,
                                         "MC event selection: require the TVX (T0 vertex) trigger"};
  Configurable<bool> mcRequireRecoINELgt0{"mcRequireRecoINELgt0", false,
                                          "MC event selection: require reconstructed INEL>0 condition"};

  // ── Histogram registry ───────────────────────────────────────────────────
  HistogramRegistry allHistograms{"allHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};

  // ============================================================
  // init()
  // ============================================================
  void init(InitContext const&)
  {
    const AxisSpec axisCentralityPercent(110, 0, 110, "FT0 centrality (%)");
    const AxisSpec axisMomentumForPID(200, 0., 10., "p (GeV/c)");
    const AxisSpec axisPtForPID(200, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisPt(numberOfPtBins, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisPseudorapidity(40, -1, 1, "#eta");
    const AxisSpec axisDCAxy(240, -0.12, 0.12, "DCA_{xy} (cm)");
    const AxisSpec axisTPCNClusters(200, 0, 200, "TPC N_{clusters}");
    const AxisSpec axisTPCNSigma(401, -10.025, 10.025, "n#sigma^{TPC}");
    const AxisSpec axisTOFNSigma(401, -10.025, 10.025, "n#sigma^{TOF}");
    const AxisSpec axisTPCdEdx(380, 10, 200, "#frac{dE}{dx}");
    const AxisSpec axisVertexZ(120, -12, 12, "v_{z} (cm)");
    const AxisSpec axisEventPlaneAngle(120, -3.14, 3.14, "#theta (rad)");
    const AxisSpec axisInvariantMass(numberOfInvMassBins, 1.4, 2.0,
                                     "M_{inv} (GeV/c^{2})");

    AxisSpec axisOccupancy = {occupancyBins, "Track occupancy [-40,100]"};
    AxisSpec axisDCAz = {dcaZMixingBins, "DCA_{z} (cm)"};

    allHistograms.add("Event/centralityVsOccupancy",
                      "Collision centrality vs track occupancy", kTH2F,
                      {axisCentralityPercent, axisOccupancy});
    allHistograms.add("Event/primaryVertexZ",
                      "Primary vertex z-position distribution", kTH1F,
                      {{100, -15., 15.}});
    if (doprocessMix || doprocessMixDF || doprocessMixepDF) {
      allHistograms.add("Event/mixingBins_centralityVsVtxZVsEventPlane",
                        "Event mixing bin occupancy: centrality vs vtxZ vs event plane", kTH3F,
                        {axisCentralityPercent, axisVertexZ, axisEventPlaneAngle});
    }

    allHistograms.add("QAbefore/trackEta", "Track pseudorapidity (before cuts)", kTH1F, {{50, -1.0, 1.0}});
    allHistograms.add("QAbefore/trackPt", "Track pT (before cuts)", kTH1F, {axisMomentumForPID});
    allHistograms.add("QAbefore/trackPhi", "Track azimuthal angle (before cuts)", kTH1F, {{72, 0, 6.2832}});
    allHistograms.add("QAbefore/trackEtaVsPhi", "Track eta vs phi (before cuts)", kTH2F, {axisPseudorapidity, {72, 0, 6.2832}});

    allHistograms.add("QAbefore/Proton/tpcNSigmaVsMomentum",
                      "TPC nSigma proton vs momentum (before PID cuts)", kTH2F,
                      {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAbefore/Proton/tofNSigmaVsMomentum",
                      "TOF nSigma proton vs momentum (before PID cuts)", kTH2F,
                      {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAbefore/Proton/tofNSigmaVsTPCNSigma",
                      "TOF nSigma vs TPC nSigma proton (before PID cuts)", kTH2F,
                      {axisTPCNSigma, axisTOFNSigma});

    allHistograms.add("QAbefore/Kaon/tpcNSigmaVsMomentum",
                      "TPC nSigma kaon vs momentum (before PID cuts)", kTH2F,
                      {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAbefore/Kaon/tofNSigmaVsMomentum",
                      "TOF nSigma kaon vs momentum (before PID cuts)", kTH2F,
                      {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAbefore/Kaon/tofNSigmaVsTPCNSigma",
                      "TOF nSigma vs TPC nSigma kaon (before PID cuts)", kTH2F,
                      {axisTPCNSigma, axisTOFNSigma});

    allHistograms.add("QAafter/Proton/ptVsCentrality",
                      "Proton pT vs centrality (after cuts)", kTH2F, {axisPtForPID, axisCentralityPercent});
    allHistograms.add("QAafter/Proton/dcaZVsPt",
                      "Proton DCAz vs pT (after cuts)", kTH2F, {axisPtForPID, axisDCAz});
    allHistograms.add("QAafter/Proton/dcaXYVsPt",
                      "Proton DCAxy vs pT (after cuts)", kTH2F, {axisPtForPID, axisDCAxy});
    allHistograms.add("QAafter/Proton/tpcDedxVsMomentum",
                      "Proton TPC dE/dx signal vs momentum (after cuts)", kTH2F,
                      {axisMomentumForPID, axisTPCdEdx});
    allHistograms.add("QAafter/Proton/tpcNSigmaVsPt",
                      "Proton TPC nSigma vs pT (after cuts)", kTH2F, {axisPtForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Proton/tpcNSigmaPionContamVsPt",
                      "Proton track: TPC nSigma pion contamination check", kTH2F,
                      {axisPtForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Proton/tpcNSigmaKaonContamVsPt",
                      "Proton track: TPC nSigma kaon contamination check", kTH2F,
                      {axisPtForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Proton/tpcNSigmaVsMomentum",
                      "Proton TPC nSigma vs total momentum (after cuts)", kTH2F,
                      {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Proton/tofNSigmaVsPt",
                      "Proton TOF nSigma vs pT (after cuts)", kTH2F, {axisPtForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Proton/tofNSigmaVsMomentum",
                      "Proton TOF nSigma vs total momentum (after cuts)", kTH2F,
                      {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Proton/tofNSigmaPionContamVsMomentum",
                      "Proton track: TOF nSigma pion contamination check", kTH2F,
                      {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Proton/tofNSigmaKaonContamVsMomentum",
                      "Proton track: TOF nSigma kaon contamination check", kTH2F,
                      {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Proton/tofNSigmaVsTPCNSigma",
                      "Proton TOF nSigma vs TPC nSigma (after cuts)", kTH2F,
                      {axisTPCNSigma, axisTOFNSigma});
    allHistograms.add("QAafter/Proton/tpcCrossedRowsVsPt",
                      "Proton TPC crossed rows vs pT", kTH2F, {axisPtForPID, {200, 0, 200}});
    allHistograms.add("QAafter/Proton/tpcClustersFoundVsPt",
                      "Proton TPC clusters found vs pT", kTH2F, {axisPtForPID, {200, 0, 200}});

    allHistograms.add("QAafter/Kaon/ptVsCentrality",
                      "Kaon pT vs centrality (after cuts)", kTH2F, {axisPtForPID, axisCentralityPercent});
    allHistograms.add("QAafter/Kaon/dcaZVsPt",
                      "Kaon DCAz vs pT (after cuts)", kTH2F, {axisPtForPID, axisDCAz});
    allHistograms.add("QAafter/Kaon/dcaXYVsPt",
                      "Kaon DCAxy vs pT (after cuts)", kTH2F, {axisPtForPID, axisDCAxy});
    allHistograms.add("QAafter/Kaon/tpcDedxVsMomentum",
                      "Kaon TPC dE/dx signal vs momentum (after cuts)", kTH2F,
                      {axisMomentumForPID, axisTPCdEdx});
    allHistograms.add("QAafter/Kaon/tpcNSigmaPionContamVsPt",
                      "Kaon track: TPC nSigma pion contamination check", kTH2F,
                      {axisPtForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Kaon/tpcNSigmaProtonContamVsMomentum",
                      "Kaon track: TPC nSigma proton contamination check", kTH2F,
                      {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Kaon/tpcNSigmaVsPt",
                      "Kaon TPC nSigma vs pT (after cuts)", kTH2F, {axisPtForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Kaon/tpcNSigmaVsMomentum",
                      "Kaon TPC nSigma vs total momentum (after cuts)", kTH2F,
                      {axisMomentumForPID, axisTPCNSigma});
    allHistograms.add("QAafter/Kaon/tofNSigmaVsPt",
                      "Kaon TOF nSigma vs pT (after cuts)", kTH2F, {axisPtForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Kaon/tofNSigmaVsMomentum",
                      "Kaon TOF nSigma vs total momentum (after cuts)", kTH2F,
                      {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Kaon/tofNSigmaPionContamVsMomentum",
                      "Kaon track: TOF nSigma pion contamination check", kTH2F,
                      {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Kaon/tofNSigmaProtonContamVsMomentum",
                      "Kaon track: TOF nSigma proton contamination check", kTH2F,
                      {axisMomentumForPID, axisTOFNSigma});
    allHistograms.add("QAafter/Kaon/tofNSigmaVsTPCNSigma",
                      "Kaon TOF nSigma vs TPC nSigma (after cuts)", kTH2F,
                      {axisTPCNSigma, axisTOFNSigma});
    allHistograms.add("QAafter/Kaon/tpcCrossedRowsVsPt",
                      "Kaon TPC crossed rows vs pT", kTH2F, {axisPtForPID, {200, 0, 200}});
    allHistograms.add("QAafter/Kaon/tpcClustersFoundVsPt",
                      "Kaon TPC clusters found vs pT", kTH2F, {axisPtForPID, {200, 0, 200}});

    if (!doprocessMC) {
      allHistograms.add("Analysis/invMass_UnlikeSign_ProtonPlusKaonMinus",
                        "Invariant mass: p^{+}K^{-} (Lambda(1520) signal)", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent, axisOccupancy});
      allHistograms.add("Analysis/invMass_UnlikeSign_ProtonMinusKaonPlus",
                        "Invariant mass: p^{-}K^{+} (anti-Lambda(1520) signal)", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent, axisOccupancy});
      allHistograms.add("Analysis/invMass_LikeSign_ProtonPlusKaonPlus",
                        "Invariant mass: p^{+}K^{+} (like-sign background)", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent, axisOccupancy});
      allHistograms.add("Analysis/invMass_LikeSign_ProtonMinusKaonMinus",
                        "Invariant mass: p^{-}K^{-} (like-sign background)", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent, axisOccupancy});
      allHistograms.add("Analysis/invMass_Rotated_ProtonPlusKaonMinus",
                        "Invariant mass: rotational background (Lambda(1520))", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent, axisOccupancy});
      allHistograms.add("Analysis/invMass_Rotated_ProtonMinusKaonPlus",
                        "Invariant mass: rotational background (anti-Lambda(1520))", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent, axisOccupancy});
      allHistograms.add("Analysis/invMass_Mixed_ProtonPlusKaonMinus",
                        "Invariant mass: mixed events p^{+}K^{-} (background)", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent, axisOccupancy});
      allHistograms.add("Analysis/invMass_Mixed_ProtonMinusKaonPlus",
                        "Invariant mass: mixed events p^{-}K^{+} (background)", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent, axisOccupancy});
      allHistograms.add("Analysis/invMass_Mixed_LikeSign_PlusPlus",
                        "Invariant mass: mixed events like-sign ++", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent, axisOccupancy});
      allHistograms.add("Analysis/invMass_Mixed_LikeSign_MinusMinus",
                        "Invariant mass: mixed events like-sign --", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent, axisOccupancy});
    }

    if (doprocessMC) {
      allHistograms.add("Event/mcEventSelectionCutflow",
                        "MC event selection cutflow (how many events pass each cut)", kTH1F, {{7, 0, 7}});
      allHistograms.add("QAChecks/protonRecoLevelPt",
                        "Reconstructed proton pT (after PID cuts) - MC", kTH1F, {axisPtForPID});
      allHistograms.add("QAChecks/kaonRecoLevelPt",
                        "Reconstructed kaon pT (after PID cuts) - MC", kTH1F, {axisPtForPID});
      allHistograms.add("QAChecks/protonGenLevelPt",
                        "Generated proton pT (truth level) - MC", kTH1F, {axisPtForPID});
      allHistograms.add("QAChecks/kaonGenLevelPt",
                        "Generated kaon pT (truth level) - MC", kTH1F, {axisPtForPID});
      allHistograms.add("Analysis/mcGenerated_Lambda1520",
                        "Generated Lambda(1520) invariant mass vs pT vs centrality", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcGenerated_AntiLambda1520",
                        "Generated anti-Lambda(1520) invariant mass vs pT vs centrality", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcReconstructed_Lambda1520",
                        "Reconstructed Lambda(1520) invariant mass vs pT vs centrality - MC", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcReconstructed_AntiLambda1520",
                        "Reconstructed anti-Lambda(1520) invariant mass vs pT vs centrality - MC", kTHnSparseF,
                        {axisInvariantMass, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcMassResolution_Lambda1520",
                        "Mass resolution (Mreco - Mgen) vs pT - Lambda(1520)", kTHnSparseF,
                        {{200, -0.05, 0.05}, axisPt, axisCentralityPercent});
      allHistograms.add("Analysis/mcMassResolution_AntiLambda1520",
                        "Mass resolution (Mreco - Mgen) vs pT - anti-Lambda(1520)", kTHnSparseF,
                        {{200, -0.05, 0.05}, axisPt, axisCentralityPercent});
    }

    if (doprocessMCGen) {
      allHistograms.add("SignalLoss/mcEventSelectionCutflow",
                        "MC event selection cutflow (signal loss study)", kTH1F, {{7, 0, 7}});
      allHistograms.add("SignalLoss/mTScaled_fromProton",
                        "mT-scaled Lambda(1520) pT from proton parent", kTHnSparseF,
                        {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromAntiProton",
                        "mT-scaled anti-Lambda(1520) pT from anti-proton parent", kTHnSparseF,
                        {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromLambda0",
                        "mT-scaled Lambda(1520) pT from Lambda0 parent", kTHnSparseF,
                        {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromAntiLambda0",
                        "mT-scaled anti-Lambda(1520) pT from anti-Lambda0 parent", kTHnSparseF,
                        {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromXiMinus",
                        "mT-scaled Lambda(1520) pT from Xi- parent", kTHnSparseF,
                        {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromXiPlus",
                        "mT-scaled anti-Lambda(1520) pT from Xi+ parent", kTHnSparseF,
                        {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromXi0",
                        "mT-scaled Lambda(1520) pT from Xi0 parent", kTHnSparseF,
                        {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromAntiXi0",
                        "mT-scaled anti-Lambda(1520) pT from anti-Xi0 parent", kTHnSparseF,
                        {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromOmegaMinus",
                        "mT-scaled Lambda(1520) pT from Omega- parent", kTHnSparseF,
                        {axisPt, axisCentralityPercent});
      allHistograms.add("SignalLoss/mTScaled_fromOmegaPlus",
                        "mT-scaled anti-Lambda(1520) pT from Omega+ parent", kTHnSparseF,
                        {axisPt, axisCentralityPercent});
    }
  }

  // ============================================================
  // passesBasicTrackSelection()
  // ============================================================
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

  // ============================================================
  // passesProtonDCASelection()
  // ============================================================
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

  // ============================================================
  // passesKaonDCASelection()
  // ============================================================
  template <typename TrackType>
  bool passesKaonDCASelection(TrackType const& track, float totalMomentum)
  {
    auto ptBinEdges = static_cast<std::vector<float>>(kaonDCAPtBinEdges);
    auto maxDCAxyValues = static_cast<std::vector<float>>(kaonMaxDCAxyPerPtBin);
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

  // ============================================================
  // passesProtonPID()
  // ============================================================
  template <typename TrackType>
  bool passesProtonPID(const TrackType& track, float totalMomentum)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};

    auto tpcMomBins = static_cast<std::vector<float>>(protonTPCPIDMomentumBins);
    auto tpcNSigCuts = static_cast<std::vector<float>>(protonTPCNSigmaCutPerBin);
    auto tofMomBins = static_cast<std::vector<float>>(protonTOFPIDMomentumBins);
    auto tofNSigCuts = static_cast<std::vector<float>>(protonTOFNSigmaCutPerBin);
    int nTPCBins = static_cast<int>(tpcMomBins.size());
    int nTOFBins = static_cast<int>(tofMomBins.size());

    float tpcNSigPion = std::abs(track.tpcNSigmaPi());
    float tpcNSigKaon = std::abs(track.tpcNSigmaKa());
    float tpcNSigProton = std::abs(track.tpcNSigmaPr());
    float tofNSigPion = std::abs(track.tofNSigmaPi());
    float tofNSigKaon = std::abs(track.tofNSigmaKa());
    float tofNSigProton = std::abs(track.tofNSigmaPr());

    float combinedNSigPion = tpcNSigPion * tpcNSigPion + tofNSigPion * tofNSigPion;
    float combinedNSigKaon = tpcNSigKaon * tpcNSigKaon + tofNSigKaon * tofNSigKaon;
    float combinedNSigProton = tpcNSigProton * tpcNSigProton + tofNSigProton * tofNSigProton;

    float circularCutSquared = combinedNSigmaCutProton * combinedNSigmaCutProton;
    float circularRejCutSquared = tofNSigmaVetoThreshold * tpcNSigmaVetoThreshold;

    if (!useTPCOnlyPID && track.hasTOF()) {
      if (track.tofNSigmaPr() < minTOFNSigmaProton)
        return false;

      if (combinedNSigmaCutProton < 0 && totalMomentum >= minTrackMomentum) {
        for (int i = 0; i < nTOFBins - 1; ++i) {
          if (totalMomentum >= tofMomBins[i] && totalMomentum < tofMomBins[i + 1] &&
              (tofNSigProton < tofNSigCuts[i] &&
               tofNSigPion > tofNSigmaVetoThreshold &&
               tofNSigKaon > tofNSigmaVetoThreshold))
            tofPIDPassed = true;
        }
        if (track.tpcNSigmaPr() < minCombinedTPCTOFNSigmaProton)
          return false;
        if (tpcNSigProton < maxTPCNSigmaProton &&
            tpcNSigPion > tpcNSigmaVetoThreshold &&
            tpcNSigKaon > tpcNSigmaVetoThreshold)
          tpcPIDPassed = true;
      }

      if ((combinedNSigmaCutProton > 0) && totalMomentum >= minTrackMomentum &&
          (combinedNSigProton < circularCutSquared &&
           combinedNSigPion > circularRejCutSquared &&
           combinedNSigKaon > circularRejCutSquared)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }

      if (totalMomentum < minTrackMomentum && tpcNSigProton < maxTPCNSigmaProton) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }

    } else {
      tofPIDPassed = true;
      if (track.tpcNSigmaPr() < minTPCNSigmaProton)
        return false;
      for (int i = 0; i < nTPCBins - 1; ++i) {
        if (totalMomentum >= tpcMomBins[i] && totalMomentum < tpcMomBins[i + 1] &&
            (tpcNSigProton < tpcNSigCuts[i] &&
             tpcNSigPion > tpcNSigmaVetoPion &&
             tpcNSigKaon > tpcNSigmaVetoKaonForProton))
          tpcPIDPassed = true;
      }
    }

    return (tpcPIDPassed && tofPIDPassed);
  }

  // ============================================================
  // passesKaonPID()
  // ============================================================
  template <typename TrackType>
  bool passesKaonPID(const TrackType& track, float totalMomentum)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};

    auto tpcMomBins = static_cast<std::vector<float>>(kaonTPCPIDMomentumBins);
    auto tpcNSigCuts = static_cast<std::vector<float>>(kaonTPCNSigmaCutPerBin);
    auto tofMomBins = static_cast<std::vector<float>>(kaonTOFPIDMomentumBins);
    auto tofNSigCuts = static_cast<std::vector<float>>(kaonTOFNSigmaCutPerBin);
    int nTPCBins = static_cast<int>(tpcMomBins.size());
    int nTOFBins = static_cast<int>(tofMomBins.size());

    float tpcNSigPion = std::abs(track.tpcNSigmaPi());
    float tpcNSigKaon = std::abs(track.tpcNSigmaKa());
    float tpcNSigProton = std::abs(track.tpcNSigmaPr());
    float tofNSigPion = std::abs(track.tofNSigmaPi());
    float tofNSigKaon = std::abs(track.tofNSigmaKa());
    float tofNSigProton = std::abs(track.tofNSigmaPr());

    float combinedNSigPion = tpcNSigPion * tpcNSigPion + tofNSigPion * tofNSigPion;
    float combinedNSigKaon = tpcNSigKaon * tpcNSigKaon + tofNSigKaon * tofNSigKaon;
    float combinedNSigProton = tpcNSigProton * tpcNSigProton + tofNSigProton * tofNSigProton;
    float circularCutSquared = combinedNSigmaCutKaon * combinedNSigmaCutKaon;
    float circularRejCutSquared = tpcNSigmaVetoOtherSpecies * tofNSigmaVetoThreshold;

    if (!useTPCOnlyPID && track.hasTOF()) {
      if (track.tofNSigmaKa() < minTOFNSigmaKaon)
        return false;

      if (combinedNSigmaCutKaon < 0 && totalMomentum >= minTrackMomentum) {
        for (int i = 0; i < nTOFBins - 1; ++i) {
          if (totalMomentum >= tofMomBins[i] && totalMomentum < tofMomBins[i + 1] &&
              (tofNSigKaon < tofNSigCuts[i] &&
               tofNSigPion > tofNSigmaVetoThreshold &&
               tofNSigProton > tofNSigmaVetoThreshold))
            tofPIDPassed = true;
        }
        if (track.tpcNSigmaKa() < minCombinedTPCTOFNSigmaKaon)
          return false;
        if (tpcNSigKaon < maxTPCNSigmaKaon &&
            tpcNSigPion > tpcNSigmaVetoThreshold &&
            tpcNSigProton > tpcNSigmaVetoThreshold)
          tpcPIDPassed = true;
      }

      if ((combinedNSigmaCutKaon > 0) && totalMomentum >= minTrackMomentum &&
          (combinedNSigKaon < circularCutSquared &&
           combinedNSigPion > circularRejCutSquared &&
           combinedNSigProton > circularRejCutSquared)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }

      if (totalMomentum < minTrackMomentum && tpcNSigKaon < maxTPCNSigmaKaon) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }

    } else {
      tofPIDPassed = true;
      if (track.tpcNSigmaKa() < minTPCNSigmaKaon)
        return false;
      for (int i = 0; i < nTPCBins - 1; ++i) {
        if (totalMomentum >= tpcMomBins[i] && totalMomentum < tpcMomBins[i + 1] &&
            (tpcNSigKaon < tpcNSigCuts[i] &&
             tpcNSigPion > tpcNSigmaVetoPionForKaon &&
             tpcNSigProton > tpcNSigmaVetoProtonForKaon))
          tpcPIDPassed = true;
      }
    }

    return (tpcPIDPassed && tofPIDPassed);
  }

  // ============================================================
  // fillInvariantMassHistograms()
  // ============================================================
  template <bool isMixedEvent, bool isMCAnalysis, typename TrackCollectionType>
  void fillInvariantMassHistograms(TrackCollectionType const& protonCandidates,
                                   TrackCollectionType const& kaonCandidates,
                                   float centralityPercent,
                                   int occupancyValue = 100)
  {
    float protonTotalMomentum = 0., kaonTotalMomentum = 0.;

    for (auto const& [protonTrack, kaonTrack] :
         soa::combinations(soa::CombinationsFullIndexPolicy(protonCandidates, kaonCandidates))) {

      if (protonTrack.index() == kaonTrack.index())
        continue;

      if (!passesBasicTrackSelection(protonTrack) || !passesBasicTrackSelection(kaonTrack))
        continue;

      auto pxProton = protonTrack.px();
      auto pyProton = protonTrack.py();
      auto pzProton = protonTrack.pz();
      auto pxKaon = kaonTrack.px();
      auto pyKaon = kaonTrack.py();
      auto pzKaon = kaonTrack.pz();

      // FIX root/entity: replaced TMath::Sqrt with RecoDecay::p
      protonTotalMomentum = RecoDecay::p(pxProton, pyProton, pzProton);
      kaonTotalMomentum = RecoDecay::p(pxKaon, pyKaon, pzKaon);

      if (!isMixedEvent) {
        auto tpcNSigProton = protonTrack.tpcNSigmaPr();
        allHistograms.fill(HIST("QAbefore/Proton/tpcNSigmaVsMomentum"), protonTotalMomentum, tpcNSigProton);
        if (protonTrack.hasTOF()) {
          auto tofNSigProton = protonTrack.tofNSigmaPr();
          allHistograms.fill(HIST("QAbefore/Proton/tofNSigmaVsMomentum"), protonTotalMomentum, tofNSigProton);
          allHistograms.fill(HIST("QAbefore/Proton/tofNSigmaVsTPCNSigma"), tpcNSigProton, tofNSigProton);
        }
        auto tpcNSigKaon = kaonTrack.tpcNSigmaKa();
        allHistograms.fill(HIST("QAbefore/Kaon/tpcNSigmaVsMomentum"), kaonTotalMomentum, tpcNSigKaon);
        if (kaonTrack.hasTOF()) {
          auto tofNSigKaon = kaonTrack.tofNSigmaKa();
          allHistograms.fill(HIST("QAbefore/Kaon/tofNSigmaVsMomentum"), kaonTotalMomentum, tofNSigKaon);
          allHistograms.fill(HIST("QAbefore/Kaon/tofNSigmaVsTPCNSigma"), tpcNSigKaon, tofNSigKaon);
        }
      }

      if (requireTOFForProton && !protonTrack.hasTOF())
        continue;
      if (requireTOFForKaon && !kaonTrack.hasTOF())
        continue;

      if (!passesProtonPID(protonTrack, protonTotalMomentum) ||
          !passesKaonPID(kaonTrack, kaonTotalMomentum))
        continue;

      if (!passesProtonDCASelection(protonTrack, protonTotalMomentum) ||
          !passesKaonDCASelection(kaonTrack, kaonTotalMomentum))
        continue;

      // FIX root/entity: replaced TMath::ACos with std::acos
      if (applyDeepAngleCut &&
          std::acos((protonTrack.pt() * kaonTrack.pt() + pzProton * pzKaon) /
                    (protonTotalMomentum * kaonTotalMomentum)) < deepAngleCutValue)
        continue;

      if constexpr (!isMixedEvent) {
        auto ptProton = protonTrack.pt();
        auto tpcNSigProton = protonTrack.tpcNSigmaPr();

        allHistograms.fill(HIST("QAafter/Proton/ptVsCentrality"), ptProton, centralityPercent);
        allHistograms.fill(HIST("QAafter/Proton/dcaZVsPt"), ptProton, protonTrack.dcaZ());
        allHistograms.fill(HIST("QAafter/Proton/dcaXYVsPt"), ptProton, protonTrack.dcaXY());
        allHistograms.fill(HIST("QAafter/Proton/tpcDedxVsMomentum"), protonTotalMomentum, protonTrack.tpcSignal());
        allHistograms.fill(HIST("QAafter/Proton/tpcNSigmaPionContamVsPt"), protonTotalMomentum, protonTrack.tpcNSigmaPi());
        allHistograms.fill(HIST("QAafter/Proton/tpcNSigmaKaonContamVsPt"), protonTotalMomentum, protonTrack.tpcNSigmaKa());
        allHistograms.fill(HIST("QAafter/Proton/tpcNSigmaVsMomentum"), protonTotalMomentum, tpcNSigProton);
        allHistograms.fill(HIST("QAafter/Proton/tpcNSigmaVsPt"), ptProton, tpcNSigProton);
        allHistograms.fill(HIST("QAafter/Proton/tpcCrossedRowsVsPt"), ptProton, protonTrack.tpcNClsCrossedRows());
        allHistograms.fill(HIST("QAafter/Proton/tpcClustersFoundVsPt"), ptProton, protonTrack.tpcNClsFound());

        if (!useTPCOnlyPID && protonTrack.hasTOF()) {
          auto tofNSigProton = protonTrack.tofNSigmaPr();
          allHistograms.fill(HIST("QAafter/Proton/tofNSigmaVsMomentum"), protonTotalMomentum, tofNSigProton);
          allHistograms.fill(HIST("QAafter/Proton/tofNSigmaVsPt"), ptProton, tofNSigProton);
          allHistograms.fill(HIST("QAafter/Proton/tofNSigmaPionContamVsMomentum"), protonTotalMomentum, protonTrack.tofNSigmaPi());
          allHistograms.fill(HIST("QAafter/Proton/tofNSigmaKaonContamVsMomentum"), protonTotalMomentum, protonTrack.tofNSigmaKa());
          allHistograms.fill(HIST("QAafter/Proton/tofNSigmaVsTPCNSigma"), tpcNSigProton, tofNSigProton);
        }

        auto ptKaon = kaonTrack.pt();
        auto tpcNSigKaon = kaonTrack.tpcNSigmaKa();

        allHistograms.fill(HIST("QAafter/Kaon/ptVsCentrality"), ptKaon, centralityPercent);
        allHistograms.fill(HIST("QAafter/Kaon/dcaZVsPt"), ptKaon, kaonTrack.dcaZ());
        allHistograms.fill(HIST("QAafter/Kaon/dcaXYVsPt"), ptKaon, kaonTrack.dcaXY());
        allHistograms.fill(HIST("QAafter/Kaon/tpcDedxVsMomentum"), kaonTotalMomentum, kaonTrack.tpcSignal());
        allHistograms.fill(HIST("QAafter/Kaon/tpcNSigmaPionContamVsPt"), kaonTotalMomentum, kaonTrack.tpcNSigmaPi());
        allHistograms.fill(HIST("QAafter/Kaon/tpcNSigmaProtonContamVsMomentum"), kaonTotalMomentum, kaonTrack.tpcNSigmaPr());
        allHistograms.fill(HIST("QAafter/Kaon/tpcNSigmaVsMomentum"), kaonTotalMomentum, tpcNSigKaon);
        allHistograms.fill(HIST("QAafter/Kaon/tpcNSigmaVsPt"), ptKaon, tpcNSigKaon);
        allHistograms.fill(HIST("QAafter/Kaon/tpcCrossedRowsVsPt"), ptKaon, kaonTrack.tpcNClsCrossedRows());
        allHistograms.fill(HIST("QAafter/Kaon/tpcClustersFoundVsPt"), ptKaon, kaonTrack.tpcNClsFound());

        if (!useTPCOnlyPID && kaonTrack.hasTOF()) {
          auto tofNSigKaon = kaonTrack.tofNSigmaKa();
          allHistograms.fill(HIST("QAafter/Kaon/tofNSigmaVsMomentum"), kaonTotalMomentum, tofNSigKaon);
          allHistograms.fill(HIST("QAafter/Kaon/tofNSigmaVsPt"), ptKaon, tofNSigKaon);
          allHistograms.fill(HIST("QAafter/Kaon/tofNSigmaPionContamVsMomentum"), kaonTotalMomentum, kaonTrack.tofNSigmaPi());
          allHistograms.fill(HIST("QAafter/Kaon/tofNSigmaProtonContamVsMomentum"), kaonTotalMomentum, kaonTrack.tofNSigmaPr());
          allHistograms.fill(HIST("QAafter/Kaon/tofNSigmaVsTPCNSigma"), tpcNSigKaon, tofNSigKaon);
        }
      }

      if (runQualityChecksOnly)
        continue;

      std::array<float, 3> momProton = {pxProton, pyProton, pzProton};
      std::array<float, 3> momKaon = {pxKaon, pyKaon, pzKaon};
      std::array<std::array<float, 3>, 2> bothMomenta = {momProton, momKaon};

      float pairInvMass = RecoDecay::m(bothMomenta, std::array{MassProton, MassKaonCharged});
      float pairPt = RecoDecay::pt(std::array{pxProton + pxKaon, pyProton + pyKaon});
      float pairRapidity = RecoDecay::y(
        std::array{pxProton + pxKaon, pyProton + pyKaon, pzProton + pzKaon},
        pairInvMass);

      if (pairRapidity < minPairRapidity || pairRapidity > maxPairRapidity)
        continue;

      if constexpr (!isMixedEvent && !isMCAnalysis) {
        if (protonTrack.sign() * kaonTrack.sign() < 0) {
          if (protonTrack.sign() > 0)
            allHistograms.fill(HIST("Analysis/invMass_UnlikeSign_ProtonPlusKaonMinus"),
                               pairInvMass, pairPt, centralityPercent, occupancyValue);
          else
            allHistograms.fill(HIST("Analysis/invMass_UnlikeSign_ProtonMinusKaonPlus"),
                               pairInvMass, pairPt, centralityPercent, occupancyValue);

          if (enableRotationalBackground) {
            for (int iRot = 0; iRot < numberOfRotations; iRot++) {
              float rotationWindowHalfWidth = o2::constants::math::PI / rotationAngleWindow;
              float rotatedKaonPhi;
              if (numberOfRotations == 1) {
                rotatedKaonPhi = o2::constants::math::PI;
              } else {
                rotatedKaonPhi = (o2::constants::math::PI - rotationWindowHalfWidth) +
                                 iRot * (2.f * rotationWindowHalfWidth / (numberOfRotations - 1));
              }

              // FIX two-pi-add-subtract: replaced manual +/-TwoPI with RecoDecay::constrainAngle
              float newKaonPhi = RecoDecay::constrainAngle(kaonTrack.phi() + rotatedKaonPhi, 0.f);

              float pxKaonRotated = kaonTrack.pt() * std::cos(newKaonPhi);
              float pyKaonRotated = kaonTrack.pt() * std::sin(newKaonPhi);

              std::array<float, 3> momProtonRot = {pxProton, pyProton, pzProton};
              std::array<float, 3> momKaonRot = {pxKaonRotated, pyKaonRotated, pzKaon};
              std::array<std::array<float, 3>, 2> rotatedMomenta = {momProtonRot, momKaonRot};

              float rotatedPairInvMass = RecoDecay::m(rotatedMomenta, std::array{MassProton, MassKaonCharged});
              float rotatedPairPt = RecoDecay::pt(std::array{pxProton + pxKaonRotated, pyProton + pyKaonRotated});
              float rotatedPairRapidity = RecoDecay::y(
                std::array{pxProton + pxKaonRotated, pyProton + pyKaonRotated, pzProton + pzKaon},
                MassLambda1520);

              if (rotatedPairRapidity < minPairRapidity || rotatedPairRapidity > maxPairRapidity)
                continue;

              if (protonTrack.sign() * kaonTrack.sign() < 0) {
                if (protonTrack.sign() > 0)
                  allHistograms.fill(HIST("Analysis/invMass_Rotated_ProtonPlusKaonMinus"),
                                     rotatedPairInvMass, rotatedPairPt, centralityPercent, occupancyValue);
                else
                  allHistograms.fill(HIST("Analysis/invMass_Rotated_ProtonMinusKaonPlus"),
                                     rotatedPairInvMass, rotatedPairPt, centralityPercent, occupancyValue);
              }
            }
          }

        } else {
          if (protonTrack.sign() > 0)
            allHistograms.fill(HIST("Analysis/invMass_LikeSign_ProtonPlusKaonPlus"),
                               pairInvMass, pairPt, centralityPercent, occupancyValue);
          else
            allHistograms.fill(HIST("Analysis/invMass_LikeSign_ProtonMinusKaonMinus"),
                               pairInvMass, pairPt, centralityPercent, occupancyValue);
        }
      }

      if constexpr (isMixedEvent) {
        if (protonTrack.sign() * kaonTrack.sign() < 0) {
          if (protonTrack.sign() > 0)
            allHistograms.fill(HIST("Analysis/invMass_Mixed_ProtonPlusKaonMinus"),
                               pairInvMass, pairPt, centralityPercent, occupancyValue);
          else
            allHistograms.fill(HIST("Analysis/invMass_Mixed_ProtonMinusKaonPlus"),
                               pairInvMass, pairPt, centralityPercent, occupancyValue);
        } else {
          if (protonTrack.sign() > 0)
            allHistograms.fill(HIST("Analysis/invMass_Mixed_LikeSign_PlusPlus"),
                               pairInvMass, pairPt, centralityPercent, occupancyValue);
          else
            allHistograms.fill(HIST("Analysis/invMass_Mixed_LikeSign_MinusMinus"),
                               pairInvMass, pairPt, centralityPercent, occupancyValue);
        }
      }

      if constexpr (isMCAnalysis) {
        if (protonTrack.sign() * kaonTrack.sign() < 0) {
          // FIX pdg/explicit-code: replaced 2212 and 321 with named constants
          if (std::abs(protonTrack.pdgCode()) != kPdgProton || std::abs(kaonTrack.pdgCode()) != kPdgKaon)
            continue;
          if (protonTrack.motherId() != kaonTrack.motherId())
            continue;
          if (protonTrack.motherPDG() != kaonTrack.motherPDG())
            continue;
          if (protonTrack.pdgCode() == 0 || kaonTrack.pdgCode() == 0)
            continue;
          if (protonTrack.motherPDG() == -1 || kaonTrack.motherPDG() == -1)
            continue;
          if (std::abs(protonTrack.motherPDG()) != pdgCodeLambda1520)
            continue;

          float trueMassOfParent = 0.;
          for (auto const& parentParticle : *mcResonanceParentTable) {
            if (parentParticle.mcParticleId() == protonTrack.motherId()) {
              std::array<float, 3> parentMom = {parentParticle.px(),
                                                parentParticle.py(),
                                                parentParticle.pz()};
              trueMassOfParent = RecoDecay::m(parentMom, parentParticle.e());
              break;
            }
          }

          float massResidual = pairInvMass - trueMassOfParent;

          if (protonTrack.motherPDG() > 0) {
            allHistograms.fill(HIST("Analysis/mcReconstructed_Lambda1520"),
                               pairInvMass, pairPt, centralityPercent);
            allHistograms.fill(HIST("Analysis/mcMassResolution_Lambda1520"),
                               massResidual, pairPt, centralityPercent);
          } else {
            allHistograms.fill(HIST("Analysis/mcReconstructed_AntiLambda1520"),
                               pairInvMass, pairPt, centralityPercent);
            allHistograms.fill(HIST("Analysis/mcMassResolution_AntiLambda1520"),
                               massResidual, pairPt, centralityPercent);
          }
        }
      }
    }
  }

  // ── Type aliases ─────────────────────────────────────────────────────────
  using ResonanceCollisionsWithEP = soa::Join<aod::ResoCollisions, aod::ResoEvtPlCollisions>;
  using ResonanceMCCollisions = soa::Join<aod::ResoCollisions, aod::ResoMCCollisions>;
  using ResonanceTrackTable = aod::ResoTracks;

  // ============================================================
  // processData()
  // ============================================================
  void processData(ResonanceCollisionsWithEP::iterator const& collision,
                   ResonanceTrackTable const& tracks)
  {
    allHistograms.fill(HIST("Event/centralityVsOccupancy"), collision.cent(), 100);
    allHistograms.fill(HIST("Event/primaryVertexZ"), collision.posZ());
    fillInvariantMassHistograms<false, false>(tracks, tracks, collision.cent());
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processData,
                 "Process real collision data (same-event analysis)", true);

  // ============================================================
  // processMC()
  // ============================================================
  void processMC(ResonanceMCCollisions::iterator const& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& tracks,
                 aod::ResoMCParents const& mcParents)
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

    if (mcRequireRecoINELgt0 && !collision.isRecINELgt0())
      return;
    allHistograms.fill(HIST("Event/mcEventSelectionCutflow"), 5);

    if (mcRequireAfterAllCuts && !collision.isInAfterAllCuts())
      return;
    allHistograms.fill(HIST("Event/mcEventSelectionCutflow"), 6);

    auto centralityPercent = collision.cent();
    allHistograms.fill(HIST("Event/centralityVsOccupancy"), centralityPercent, 100);
    allHistograms.fill(HIST("Event/primaryVertexZ"), collision.posZ());

    mcResonanceParentTable = &mcParents;
    fillInvariantMassHistograms<false, true>(tracks, tracks, centralityPercent);

    // FIX const-ref-in-for-loop: added const to range-based for loop
    for (const auto& track : tracks) {
      allHistograms.fill(HIST("QAbefore/trackEta"), track.eta());
      allHistograms.fill(HIST("QAbefore/trackPt"), track.pt());
      allHistograms.fill(HIST("QAbefore/trackPhi"), track.phi());
      allHistograms.fill(HIST("QAbefore/trackEtaVsPhi"), track.eta(), track.phi());

      // FIX pdg/explicit-code: replaced 321 and 2212 with named constants
      if (std::abs(track.pdgCode()) == kPdgKaon)
        allHistograms.fill(HIST("QAChecks/kaonGenLevelPt"), track.pt());
      if (std::abs(track.pdgCode()) == kPdgProton)
        allHistograms.fill(HIST("QAChecks/protonGenLevelPt"), track.pt());

      if (!passesBasicTrackSelection(track))
        continue;

      // FIX root/entity: replaced TMath::Sqrt with RecoDecay::p
      float totalMom = RecoDecay::p(track.px(), track.py(), track.pz());

      if (passesKaonPID(track, totalMom) && std::abs(track.pdgCode()) == kPdgKaon)
        allHistograms.fill(HIST("QAChecks/kaonRecoLevelPt"), track.pt());
      if (passesProtonPID(track, totalMom) && std::abs(track.pdgCode()) == kPdgProton)
        allHistograms.fill(HIST("QAChecks/protonRecoLevelPt"), track.pt());
    }

    // FIX const-ref-in-for-loop: added const to range-based for loop
    for (const auto& parentParticle : mcParents) {
      if (std::abs(parentParticle.pdgCode()) != pdgCodeLambda1520)
        continue;
      if (parentParticle.y() < minPairRapidity || parentParticle.y() > maxPairRapidity)
        continue;

      // FIX pdg/explicit-code: replaced 2212 and 321 with named constants
      bool hasProtonDaughter = (std::abs(parentParticle.daughterPDG1()) == kPdgProton ||
                                std::abs(parentParticle.daughterPDG2()) == kPdgProton);
      bool hasKaonDaughter = (std::abs(parentParticle.daughterPDG1()) == kPdgKaon ||
                              std::abs(parentParticle.daughterPDG2()) == kPdgKaon);
      if (!hasProtonDaughter || !hasKaonDaughter)
        continue;

      std::array<float, 3> parentMom = {parentParticle.px(), parentParticle.py(), parentParticle.pz()};
      float parentMass = RecoDecay::m(parentMom, parentParticle.e());

      if (parentParticle.pdgCode() > 0)
        allHistograms.fill(HIST("Analysis/mcGenerated_Lambda1520"),
                           parentMass, parentParticle.pt(), centralityPercent);
      else
        allHistograms.fill(HIST("Analysis/mcGenerated_AntiLambda1520"),
                           parentMass, parentParticle.pt(), centralityPercent);
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processMC,
                 "Process Monte Carlo simulated events", false);

  // ============================================================
  // processMCGen()
  // ============================================================
  void processMCGen(ResonanceMCCollisions::iterator const& collision,
                    aod::ResoMCParents const& mcParents)
  {
    float centralityPercent = collision.cent();

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
    if (mcRequireRecoINELgt0 && !collision.isRecINELgt0())
      return;
    allHistograms.fill(HIST("SignalLoss/mcEventSelectionCutflow"), 5);
    if (mcRequireAfterAllCuts && !collision.isInAfterAllCuts())
      return;
    allHistograms.fill(HIST("SignalLoss/mcEventSelectionCutflow"), 6);

    // FIX const-ref-in-for-loop: added const to range-based for loop
    for (const auto& parentParticle : mcParents) {
      if (parentParticle.y() < minPairRapidity || parentParticle.y() > maxPairRapidity)
        continue;

      int pdgCode = parentParticle.pdgCode();
      float parentPt = parentParticle.pt();
      double mTScaledPtSquared = -1.0;

      std::array<float, 3> parentMom = {parentParticle.px(), parentParticle.py(), parentParticle.pz()};
      float parentMass = RecoDecay::m(parentMom, parentParticle.e());

      auto computeMtScaledPtSquared = [&]() -> double {
        return (parentPt * parentPt) + (parentMass * parentMass) -
               (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
      };

      // FIX pdg/explicit-code: replaced all bare PDG integers with named constants
      if (pdgCode == kPdgProton) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromProton"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == -kPdgProton) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromAntiProton"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == kPdgLambda0) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromLambda0"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == -kPdgLambda0) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromAntiLambda0"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == kPdgXiMinus) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromXiMinus"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == -kPdgXiMinus) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromXiPlus"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == kPdgXi0) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromXi0"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == -kPdgXi0) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromAntiXi0"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == kPdgOmegaMinus) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromOmegaMinus"), std::sqrt(mTScaledPtSquared), centralityPercent);
      } else if (pdgCode == -kPdgOmegaMinus) {
        mTScaledPtSquared = computeMtScaledPtSquared();
        if (mTScaledPtSquared > 0)
          allHistograms.fill(HIST("SignalLoss/mTScaled_fromOmegaPlus"), std::sqrt(mTScaledPtSquared), centralityPercent);
      }
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processMCGen,
                 "Generator-level MC signal loss study (mT scaling)", false);

  // ── Event-mixing binning types ───────────────────────────────────────────
  using MixingBinningVtxZAndCentrality =
    ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;

  // ============================================================
  // processMix()
  // FIX const-ref-in-process: collisions is now const&
  // ============================================================
  void processMix(ResonanceCollisionsWithEP const& collisions,
                  ResonanceTrackTable const& tracks)
  {
    LOGF(debug, "Event mixing started");
    MixingBinningVtxZAndCentrality mixingBins{{vertexZMixingBins, centralityMixingBins}, true};
    auto trackPool = std::make_tuple(tracks);
    SameKindPair<ResonanceCollisionsWithEP, ResonanceTrackTable,
                 MixingBinningVtxZAndCentrality>
      eventPairs{mixingBins, numberOfEventsToMix, -1, collisions, trackPool, &sliceCache};

    // FIX const-ref-in-for-loop: added const to structured binding
    for (const auto& [col1, tracks1, col2, tracks2] : eventPairs) {
      allHistograms.fill(HIST("Event/mixingBins_centralityVsVtxZVsEventPlane"),
                         col1.cent(), col1.posZ(), col1.evtPl());
      fillInvariantMassHistograms<true, false>(tracks1, tracks2, col1.cent());
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processMix,
                 "Event mixing for background estimation (standard format)", true);

  // ── Merged-DF type aliases ───────────────────────────────────────────────
  Preslice<aod::ResoTrackDFs> tracksPerMergedDFCollision = aod::resodaughter::resoCollisionDFId;
  using MergedDFCollisions = aod::ResoCollisionDFs;
  using MergedDFTracks = aod::ResoTrackDFs;

  // ============================================================
  // processDatadf()
  // ============================================================
  void processDatadf(MergedDFCollisions::iterator const& collision,
                     MergedDFTracks const& tracks)
  {
    if (doprocessData)
      LOG(error) << "Disable processData() first when using processDatadf()!";

    auto occupancyValue = 100;
    if (applyOccupancyInTimeRangeCut)
      occupancyValue = collision.trackOccupancyInTimeRange();

    allHistograms.fill(HIST("Event/centralityVsOccupancy"), collision.cent(), occupancyValue);
    fillInvariantMassHistograms<false, false>(tracks, tracks, collision.cent(), occupancyValue);
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processDatadf,
                 "Process real data in merged derived-data (DF) format", false);

  // ============================================================
  // processMixDF()
  // FIX const-ref-in-process: collisions is now const&
  // ============================================================
  using MixingBinningDF = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;

  void processMixDF(MergedDFCollisions const& collisions, MergedDFTracks const& tracks)
  {
    if (doprocessMix)
      LOG(fatal) << "Disable processMix() first when using processMixDF()!";
    LOGF(debug, "Event mixing (DF format) started");

    MixingBinningDF mixingBins{{vertexZMixingBins, centralityMixingBins}, true};
    auto trackPool = std::make_tuple(tracks);
    SameKindPair<MergedDFCollisions, MergedDFTracks, MixingBinningDF>
      eventPairs{mixingBins, numberOfEventsToMix, -1, collisions, trackPool, &sliceCache};

    // FIX const-ref-in-for-loop: added const to structured binding
    for (const auto& [col1, tracks1, col2, tracks2] : eventPairs) {
      auto occupancyValue = 100;
      if (applyOccupancyInTimeRangeCut)
        occupancyValue = col1.trackOccupancyInTimeRange();

      allHistograms.fill(HIST("Event/mixingBins_centralityVsVtxZVsEventPlane"),
                         col1.cent(), col1.posZ(), col1.evtPl());
      fillInvariantMassHistograms<true, false>(tracks1, tracks2, col1.cent(), occupancyValue);
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processMixDF,
                 "Event mixing for DF-format data", false);

  // ============================================================
  // processMixepDF()
  // FIX const-ref-in-process: collisions is now const&
  // ============================================================
  using MixingBinningWithEventPlane =
    ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent, aod::resocollision::EvtPl>;

  void processMixepDF(MergedDFCollisions const& collisions, MergedDFTracks const& tracks)
  {
    if (doprocessMix || doprocessMixDF)
      LOG(fatal) << "Disable processMix() or processMixDF() first!";
    LOGF(debug, "Event-plane-dependent event mixing (DF format) started");

    MixingBinningWithEventPlane mixingBins{
      {vertexZMixingBins, centralityMixingBins, eventPlaneMixingBins}, true};
    auto trackPool = std::make_tuple(tracks);
    SameKindPair<MergedDFCollisions, MergedDFTracks, MixingBinningWithEventPlane>
      eventPairs{mixingBins, numberOfEventsToMix, -1, collisions, trackPool, &sliceCache};

    // FIX const-ref-in-for-loop: added const to structured binding
    for (const auto& [col1, tracks1, col2, tracks2] : eventPairs) {
      allHistograms.fill(HIST("Event/mixingBins_centralityVsVtxZVsEventPlane"),
                         col1.cent(), col1.posZ(), col1.evtPl());
      fillInvariantMassHistograms<true, false>(tracks1, tracks2, col1.cent());
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processMixepDF,
                 "Event-plane dependent event mixing for DF-format data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Lambda1520Analysispo>(cfgc)};
}
