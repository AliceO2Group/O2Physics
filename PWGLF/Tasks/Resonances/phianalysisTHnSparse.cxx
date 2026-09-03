// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
/// \file phianalysisTHnSparse.cxx
/// \brief Analysis of phi resonance using THnSparse histograms.
/// \author Veronika Barbasova (veronika.barbasova@cern.ch)

#include "PWGLF/Utils/rsnOutput.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/GenVector/AxisAngle.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TString.h>

#include <fmt/format.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PhianalysisTHnSparse {

  SliceCache cache;

  struct : ConfigurableGroup {
    Configurable<bool> produceMC{"produceMC", false, "Produce True and Gen histograms."};
    Configurable<bool> produceLikesign{"produceLikesign", false, "Produce Like sign histograms."};
    Configurable<std::string> eventMixing{"eventMixing", "none", "Produce Event Mixing histograms of type."};
    Configurable<bool> produceRotational{"produceRotational", false, "Produce Rotational histograms."};
  } produce;

  Configurable<int> daughterPos{"daughterPos", 3, "Particle type of the positive daughter according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> daughterNeg{"daughterNeg", 3, "Particle type of the negative daughter according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> motherPDG{"motherPDG", 333, "PDG code of mother particle."};
  Configurable<int> daughterPosPDG{"daughterPosPDG", 321, "PDG code of positive daughter particle."};
  Configurable<int> daughterNegPDG{"daughterNegPDG", 321, "PDG code of negative daughter particle."};

  struct : ConfigurableGroup {
    Configurable<bool> isTriggerTVX{"isTriggerTVX", false, "Apply IsTriggerTVX cut."};
    Configurable<bool> noTimeFrameBorder{"noTimeFrameBorder", false, "Apply NoTimeFrameBorder cut."};
    Configurable<bool> noITSROFrameBorder{"noITSROFrameBorder", false, "Apply NoITSROFrameBorder cut."};
    Configurable<bool> sel8{"sel8", false, "Apply Sel8 cut."};
    Configurable<bool> inelGt0{"inelGt0", false, "Select events with INEL>0."};
    Configurable<float> vzCut{"vzCut", 10.0f, "Cut: Maximal value of Z vertex position."};
    Configurable<bool> noSameBunchPileup{"noSameBunchPileup", false, "Apply no same bunch pileup cut."};
    Configurable<bool> isVertexITSTPC{"isVertexITSTPC", false, "Apply IsVertexITSTPC cut."};
    Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", false, "Apply IsGoodZvtxFT0vsPV cut."};
  } eventCuts;

  struct : ConfigurableGroup {
    Configurable<float> pt{"pt", 0.15f, "Cut: Minimal value of tracks pt."};
    Configurable<float> etatrack{"etatrack", 1.0f, "Cut: Maximal value of tracks eta."};
    Configurable<float> dcaXY{"dcaXY", 1.0f, "Cut: Maximal value of tracks DCA XY."};
    Configurable<float> dcaZ{"dcaZ", 1.0f, "Cut: Maximal value of tracks DCA Z."};
    Configurable<float> tpcnSigmaPos{"tpcnSigmaPos", 10.0f, "Cut: Maximal value of TPC NSigma of the positive particle."};
    Configurable<float> tpcnSigmaNeg{"tpcnSigmaNeg", 10.0f, "Cut: Maximal value of TPC NSigma of the negative particle."};
    Configurable<bool> tpcPidOnly{"tpcPidOnly", false, "Use TPC only for PID."};
    Configurable<float> combinedNSigma{"combinedNSigma", 3.0f, "Cut: Maximal value of NSigma for combined TPC and TOF NSigma cut."};
    Configurable<float> ptTOFThreshold{"ptTOFThreshold", 0.5f, "Cut: Minimal value of tracks pt for using TOF PID."};
    Configurable<int> tpcNClsFound{"tpcNClsFound", 155, "Cut: Minimal value of found TPC clusters"};
    Configurable<int> tpcNClsCrossedRows{"tpcNClsCrossedRows", 155, "Cut: Minimal value of crossed rows in TPC"};
    Configurable<bool> globalTrack{"globalTrack", false, "Use isGlobalTrack track selection."};
    Configurable<bool> primaryTrack{"primaryTrack", false, "Use isPrimaryTrack track selection."};
    Configurable<bool> pvContributor{"pvContributor", false, "Use isPVContributor track selection."};
    Configurable<float> rapidity{"rapidity", 0.5f, "Cut: Maximal value of particle rapidity."};
  } trackCuts;

  Configurable<std::vector<std::string>> sparseAxes{"sparseAxes", std::vector<std::string>{o2::analysis::rsn::pair_axis::names}, "Axes."};
  Configurable<std::vector<std::string>> sysAxes{"sysAxes", std::vector<std::string>{o2::analysis::rsn::systematic_axis::names}, "Axes."};

  ConfigurableAxis invaxis{"invaxis", {130, 0.97, 1.1}, "Invariant mass axis binning."};
  ConfigurableAxis ptaxis{"ptaxis", {20, 0., 20.}, "Pt axis binning."};
  ConfigurableAxis vzaxis{"vzaxis", {40, -20., 20.}, "Z vertex position axis binning."};
  ConfigurableAxis multiplicityaxis{"multiplicityaxis", {50, 0., 5000.}, "Multiplicity axis binning."};
  ConfigurableAxis centralityaxis{"centralityaxis", {20, 0., 100.}, "Centrality axis binning."};
  ConfigurableAxis etaaxis{"etaaxis", {16., -1.0 * static_cast<float>(trackCuts.etatrack), static_cast<float>(trackCuts.etatrack)}, "Pseudorapidity axis binning."};
  ConfigurableAxis rapidityaxis{"rapidityaxis", {10., -1.0 * static_cast<float>(trackCuts.rapidity), static_cast<float>(trackCuts.rapidity)}, "Rapidity axis binning."};
  ConfigurableAxis nsigmaaxisPos{"nsigmaaxisPos", {1, -static_cast<float>(trackCuts.tpcnSigmaPos), static_cast<float>(trackCuts.tpcnSigmaPos)}, "NSigma of positive particle axis binning in THnSparse."};
  ConfigurableAxis nsigmaaxisNeg{"nsigmaaxisNeg", {1, -static_cast<float>(trackCuts.tpcnSigmaNeg), static_cast<float>(trackCuts.tpcnSigmaNeg)}, "NSigma of negative particle axis binning in THnSparse."};

  // mixing
  using BinningTypeVzMu = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
  using BinningTypeVzCe = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  Configurable<int> nMixedEvents{"nMixedEvents", 5, "Number of events that should be mixed."};
  ConfigurableAxis axisVertexMixing{"axisVertexMixing", {5, -10, 10}, "Z vertex axis binning for mixing"};
  ConfigurableAxis axisMultiplicityMixing{"axisMultiplicityMixing", {5, 0, 5000}, "FT0M amplitude binning for event mixing."};
  ConfigurableAxis axisCentralityMixing{"axisCentralityMixing", {10, 0, 100}, "FT0M centrality percentile binning for event mixing."};

  // rotational
  Configurable<int> nRotations{"nRotations", 1, "Number of rotations for rotational background estimation."};
  Configurable<int> startingAngle{"startingAngle", 0, "Starting angle for rotational background estimation."};

  // other axes
  ConfigurableAxis axisNch{"axisNch", {1000, 0.0, +1000.0}, "Number of charged particles."};
  ConfigurableAxis axisResolutionPt{"axisResolutionPt", {1001, -1.0, +1.0}, "Resolution of Pt."};
  ConfigurableAxis axisResolutionPtPhi{"axisResolutionPtPhi", {1001, -0.01, +0.01}, "Resolution of Pt and Phi."};
  ConfigurableAxis axisResolutionMass{"axisResolutionMass", {1001, -0.01, +0.01}, "Resolution of Mass."};
  ConfigurableAxis axisResolutionVz{"axisResolutionVz", {1001, -3.0, +3.0}, "Resolution of Vz."};
  ConfigurableAxis axisQAPt{"axisQAPt", {15, 0.0, 15.0}, "QA Pt axis binning."};
  ConfigurableAxis axisQAMult{"axisQAMult", {10, 0.0, 100.0}, "QA Multiplicity axis binning."};
  ConfigurableAxis axisQACent{"axisQACent", {101, 0.0f, 101.0f}, "QA Centrality axis binning."};

  // Axes specifications
  AxisSpec vzQAaxis = {200, -20., 20., "V_{z} (cm)"};
  AxisSpec dcaXYQAaxis = {200, -0.5, 0.5, "DCA_{xy} (cm)"};
  AxisSpec dcaZQAaxis = {200, -0.5, 0.5, "DCA_{z} (cm)"};
  AxisSpec etaQAaxis = {200, -1.0, 1.0, "#eta"};
  AxisSpec rapidityQAaxis = {200, -1.0, 1.0, "y"};
  AxisSpec tpcNClsQAaxis = {200, 0., 200., "TPC NClusters"};
  AxisSpec nSigmaTPCQAaxis = {200, -10., 10., "n#sigma_{TPC} K^{#pm}"};
  AxisSpec nSigmaTOFQAaxis = {200, -10., 10., "n#sigma_{TOF} K^{#pm}"};
  AxisSpec pQAaxis = {1490, 0.1, 15.0, "p (GeV/c)"};
  AxisSpec dEdxQAaxis = {2000, 0., 200., "dE/dx (a.u.)"};
  AxisSpec betaQAaxis = {700, 0.5, 1.2, "#beta"};
  AxisSpec dPhiQAaxis = {100, -o2::constants::math::TwoPI, o2::constants::math::TwoPI, "#Delta#phi (rad)"};
  AxisSpec dThetaQAaxis = {100, -o2::constants::math::PI, o2::constants::math::PI, "#Delta#theta (rad)"};
  AxisSpec dEtaQAaxis = {200, -1.0, 1.0, "#Delta#eta"};

  HistogramRegistry registry{"registry"};
  o2::analysis::rsn::Output* rsnOutput = nullptr;

  Service<o2::framework::O2DatabasePDG> pdg{};

  float massPos = o2::track::PID::getMass(3);
  float massNeg = o2::track::PID::getMass(3);
  int pion = 2;
  int kaon = 3;
  int proton = 4;
  double* pointPair = nullptr;
  ROOT::Math::PxPyPzMVector d1, d2, mother, motherGen;
  bool dataQA = false;
  int dauSize = 2;
  rsn::MixingType mixingType = rsn::MixingType::none;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Ms>;
  using EventCandidate = EventCandidates::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTOFbeta>;

  using EventCandidatesMC = soa::Join<EventCandidates, aod::McCollisionLabels>;
  using TrackCandidatesMC = soa::Join<TrackCandidates, aod::McTrackLabels>;

  using EventCandidatesMCGen = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Ms, aod::McCollisionLabels>;
  using McCollisionMults = soa::Join<aod::McCollisions, aod::MultMCExtras>;
  using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  Partition<TrackCandidates> positive = (aod::track::signed1Pt > 0.0f);
  Partition<TrackCandidates> negative = (aod::track::signed1Pt < 0.0f);

  Partition<TrackCandidatesMC> positiveMC = (aod::track::signed1Pt > 0.0f);
  Partition<TrackCandidatesMC> negativeMC = (aod::track::signed1Pt < 0.0f);

  void init(o2::framework::InitContext&)
  {
    // defined in DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h
    massPos = o2::track::PID::getMass(static_cast<int>(daughterPos));
    massNeg = o2::track::PID::getMass(static_cast<int>(daughterNeg));
    LOGF(info, "Initializing particle masses: ");
    LOGF(info, "  Positive: %d, mass: %f", static_cast<int>(daughterPos), massPos);
    LOGF(info, "  Negative: %d, mass: %f", static_cast<int>(daughterNeg), massNeg);

    AxisSpec centQAAxis = {axisQACent, "FT0M (%)"};
    AxisSpec nchQAAxis = {axisNch, "N_{ch}"};
    AxisSpec multQAAxis = {axisQAMult, "FT0M (%)"};
    AxisSpec ptQAAxis = {axisQAPt, "p_{T} (GeV/c)"};

    // Sparse axes
    AxisSpec invAxis = {invaxis, "Inv. mass (GeV/c^{2})", "im"};
    AxisSpec ptAxis = {ptaxis, "p_{T} (GeV/c)", "pt"};
    AxisSpec muAxis = {multiplicityaxis, "FT0M (Ampl.)", "mu"};
    AxisSpec mumAxis = {multiplicityaxis, "FT0M (Ampl.)", "mum"};
    AxisSpec ceAxis = {centralityaxis, "FT0M (%)", "ce"};
    AxisSpec cemAxis = {centralityaxis, "FT0M (%)", "cem"};
    AxisSpec etaAxis = {etaaxis, "#eta", "eta"};
    AxisSpec yAxis = {rapidityaxis, "y", "y"};
    AxisSpec nsAxisPos = {nsigmaaxisPos, fmt::format("nSigma of positive particle ({})", massPos), "ns1"};
    AxisSpec nsAxisNeg = {nsigmaaxisNeg, fmt::format("nSigma of negative particle ({})", massNeg), "ns2"};
    AxisSpec vzAxis = {vzaxis, "V_{z} (cm)", "vz"};
    AxisSpec vzmAxis = {vzaxis, "V_{z} (cm)", "vzm"};

    // Systematics axes
    std::vector<double> tpcNClsFoundBins = {65., 66., 67., 68., 69., 70., 71., 72., 73.};
    AxisSpec tpcNClsFoundAxis = {tpcNClsFoundBins, "TPC NCl", "ncl"};

    // All axes has to have same order as defined enum o2::analysis::rsn::PairAxisType (name from AxisSpec is taken to compare in o2::analysis::rsn::Output::init())
    std::vector<AxisSpec> allAxes = {invAxis, ptAxis, muAxis, ceAxis, nsAxisPos, nsAxisNeg, etaAxis, yAxis, vzAxis, mumAxis, cemAxis, vzmAxis};
    std::vector<AxisSpec> allAxesSys = {tpcNClsFoundAxis};

    mixingType = rsn::mixingTypeName(static_cast<std::string>(produce.eventMixing));

    pointPair = new double[static_cast<int>(o2::analysis::rsn::PairAxisType::unknown)];
    rsnOutput = new o2::analysis::rsn::OutputSparse();
    rsnOutput->init(sparseAxes, allAxes, sysAxes, allAxesSys, static_cast<bool>(produce.produceMC), mixingType, static_cast<bool>(produce.produceLikesign), static_cast<bool>(produce.produceRotational), &registry);

    // Print summary of configuration
    LOGF(info, "=== PhianalysisTHnSparse configuration summary ===");
    LOGF(info, "produceMC: %s", static_cast<bool>(produce.produceMC) ? "true" : "false");
    LOGF(info, "produceLikesign: %s", static_cast<bool>(produce.produceLikesign) ? "true" : "false");
    LOGF(info, "produceRotational: %s", static_cast<bool>(produce.produceRotational) ? "true" : "false");
    LOGF(info, "eventMixing: %s", static_cast<std::string>(produce.eventMixing).c_str());
    LOGF(info, "inelGt0: %s", static_cast<bool>(eventCuts.inelGt0) ? "true" : "false");
    LOGF(info, "noSameBunchPileup: %s", static_cast<bool>(eventCuts.noSameBunchPileup) ? "true" : "false");
    LOGF(info, "isVertexITSTPC: %s", static_cast<bool>(eventCuts.isVertexITSTPC) ? "true" : "false");
    LOGF(info, "isGoodZvtxFT0vsPV: %s", static_cast<bool>(eventCuts.isGoodZvtxFT0vsPV) ? "true" : "false");
    LOGF(info, "vzCut: %.2f", static_cast<float>(eventCuts.vzCut));
    LOGF(info, "daughterPos: %d (PDG: %d)", static_cast<int>(daughterPos), static_cast<int>(daughterPosPDG));
    LOGF(info, "daughterNeg: %d (PDG: %d)", static_cast<int>(daughterNeg), static_cast<int>(daughterNegPDG));
    LOGF(info, "motherPDG: %d", static_cast<int>(motherPDG));
    LOGF(info, "pt (min): %.2f", static_cast<float>(trackCuts.pt));
    LOGF(info, "eta (max): %.2f", static_cast<float>(trackCuts.etatrack));
    LOGF(info, "dcaXY: %.2f", static_cast<float>(trackCuts.dcaXY));
    LOGF(info, "dcaZ: %.2f", static_cast<float>(trackCuts.dcaZ));
    LOGF(info, "tpcnSigmaPos: %.2f", static_cast<float>(trackCuts.tpcnSigmaPos));
    LOGF(info, "tpcnSigmaNeg: %.2f", static_cast<float>(trackCuts.tpcnSigmaNeg));
    LOGF(info, "tpcPidOnly: %s", static_cast<bool>(trackCuts.tpcPidOnly) ? "true" : "false");
    LOGF(info, "combinedNSigma: %.2f", static_cast<float>(trackCuts.combinedNSigma));
    LOGF(info, "ptTOFThreshold: %.2f", static_cast<float>(trackCuts.ptTOFThreshold));
    LOGF(info, "tpcNClsFound: %d", static_cast<int>(trackCuts.tpcNClsFound));
    LOGF(info, "tpcNClsCrossedRows: %d", static_cast<int>(trackCuts.tpcNClsCrossedRows));
    LOGF(info, "globalTrack: %s", static_cast<bool>(trackCuts.globalTrack) ? "true" : "false");
    LOGF(info, "primaryTrack: %s", static_cast<bool>(trackCuts.primaryTrack) ? "true" : "false");
    LOGF(info, "pvContributor: %s", static_cast<bool>(trackCuts.pvContributor) ? "true" : "false");
    LOGF(info, "rapidity: %.2f", static_cast<float>(trackCuts.rapidity));
    LOGF(info, "nMixedEvents: %d", static_cast<int>(nMixedEvents));
    LOGF(info, "nRotations: %d", static_cast<int>(nRotations));
    LOGF(info, "startingAngle: %d", static_cast<int>(startingAngle));
    LOGF(info, "sparseAxes: ");
    for (const auto& axis : static_cast<std::vector<std::string>>(sparseAxes)) {
      LOGF(info, "  %s", axis.c_str());
    }
    LOGF(info, "sysAxes: ");
    for (const auto& axis : static_cast<std::vector<std::string>>(sysAxes)) {
      LOGF(info, "  %s", axis.c_str());
    }
    LOGF(info, "===============================================");

    // ------------------- Event QA -------------------
    registry.add("QA/Event/hSelection", "Event selection statistics", kTH1D, {{11, 0.0f, 11.0f}});
    auto hEvent = registry.get<TH1>(HIST("QA/Event/hSelection"));
    hEvent->GetXaxis()->SetBinLabel(1, "all events");
    hEvent->GetXaxis()->SetBinLabel(2, "isTriggerTVX");
    hEvent->GetXaxis()->SetBinLabel(3, "noTimeFrameBorder");
    hEvent->GetXaxis()->SetBinLabel(4, "noITSROFrameBorder");
    hEvent->GetXaxis()->SetBinLabel(5, "sel8");
    hEvent->GetXaxis()->SetBinLabel(6, "IsVertexITSTPC");
    hEvent->GetXaxis()->SetBinLabel(7, "noSameBunchPileup");
    hEvent->GetXaxis()->SetBinLabel(8, "IsGoodZvtxFT0vsPV");
    hEvent->GetXaxis()->SetBinLabel(9, Form("|V_{z}| < %0.0f cm", static_cast<float>(eventCuts.vzCut)));
    hEvent->GetXaxis()->SetBinLabel(10, "INEL");
    hEvent->GetXaxis()->SetBinLabel(11, "INEL>0");
    hEvent->SetMinimum(0.1);

    registry.add("QA/Event/hVtxZ", "Vertex position along the z-axis", kTH1F, {vzQAaxis});
    auto hVtxZ = registry.get<TH1>(HIST("QA/Event/hVtxZ"));

    registry.add("QA/Event/hCent", "FT0M (%)", kTH1F, {{101, 0., 101.}});
    auto hCent = registry.get<TH1>(HIST("QA/Event/hCent"));
    hCent->GetXaxis()->SetTitle("FT0M (%)");

    registry.add("QA/Event/hMult", "Amplitude of non-zero channels in the FT0A + FT0C) ", kTH1F, {{300, 0., 30000.}});
    auto hMult = registry.get<TH1>(HIST("QA/Event/hMult"));
    hMult->GetXaxis()->SetTitle("FT0M Ampl.");

    registry.add("QA/Event/hCentNch", "Event centrality vs multiplicity", kTH2F, {centQAAxis, nchQAAxis});

    // ----------------------- Track QA -----------------------
    registry.add("QA/Track/hSelection", "Track selection statistics", kTH1D, {{9, 0.0f, 9.0f}});
    auto hTrack = registry.get<TH1>(HIST("QA/Track/hSelection"));
    hTrack->GetXaxis()->SetBinLabel(1, "all tracks");
    hTrack->GetXaxis()->SetBinLabel(2, Form("pT > %.2f", static_cast<float>(trackCuts.pt)));
    hTrack->GetXaxis()->SetBinLabel(3, Form("eta < %.1f", static_cast<float>(trackCuts.etatrack)));
    hTrack->GetXaxis()->SetBinLabel(4, "DCA cuts");
    hTrack->GetXaxis()->SetBinLabel(5, "PID cuts");
    hTrack->GetXaxis()->SetBinLabel(6, Form("tpcNClsFound > %d", static_cast<int>(trackCuts.tpcNClsFound)));
    hTrack->GetXaxis()->SetBinLabel(7, Form("tpcNClsCrossedRows > %d", static_cast<int>(trackCuts.tpcNClsCrossedRows)));
    hTrack->GetXaxis()->SetBinLabel(8, Form("%s", static_cast<bool>(trackCuts.globalTrack) ? "isGlobalTrack" : "isPrimaryTrack"));
    hTrack->GetXaxis()->SetBinLabel(9, "isPVContributor");
    hTrack->SetMinimum(0.1);

    registry.add("QA/Track/hRapidity", "Rapidity distribution of Tracks", kTH3F, {ptQAAxis, multQAAxis, rapidityQAaxis});
    registry.add("QA/Track/hEta", "Pseudorapidity distribution of Tracks", kTH3F, {ptQAAxis, multQAAxis, etaQAaxis});
    registry.add("QA/Track/hTPCNClsFound", "Number of found TPC clusters of Tracks", kTH3F, {ptQAAxis, multQAAxis, tpcNClsQAaxis});
    registry.add("QA/Track/hTPCNClsCrossedRows", "Number of crossed rows in TPC of Tracks", kTH3F, {ptQAAxis, multQAAxis, tpcNClsQAaxis});
    registry.add("QA/Track/hDCAxy", "Distribution of DCA_{xy} of Tracks", kTH3F, {ptQAAxis, multQAAxis, dcaXYQAaxis});
    registry.add("QA/Track/hDCAz", "Distribution of DCA_{z} of Tracks", kTH3F, {ptQAAxis, multQAAxis, dcaZQAaxis});
    registry.add("QA/Track/hPt", "Distribution of p_{T} of Tracks", kTH2F, {ptQAAxis, multQAAxis});

    registry.add("QA/Kaon/hRapidity", "Rapidity distribution of K^{+} and K^{-}", kTH3F, {ptQAAxis, multQAAxis, rapidityQAaxis});
    registry.add("QA/Kaon/hEta", "Pseudorapidity distribution of K^{+} and K^{-}", kTH3F, {ptQAAxis, multQAAxis, etaQAaxis});
    registry.add("QA/Kaon/hTPCNClsFound", "Number of found TPC clusters of K^{+} and K^{-}", kTH3F, {ptQAAxis, multQAAxis, tpcNClsQAaxis});
    registry.add("QA/Kaon/hTPCNClsCrossedRows", "Number of crossed rows in TPC of K^{+} and K^{-}", kTH3F, {ptQAAxis, multQAAxis, tpcNClsQAaxis});
    registry.add("QA/Kaon/hDCAxy", "Distribution of DCA_{xy} of K^{+} and K^{-}", kTH3F, {ptQAAxis, multQAAxis, dcaXYQAaxis});
    registry.add("QA/Kaon/hDCAz", "Distribution of DCA_{z} of K^{+} and K^{-}", kTH3F, {ptQAAxis, multQAAxis, dcaZQAaxis});
    registry.add("QA/Kaon/hPt", "Distribution of p_{T} of K^{+} and K^{-}", kTH2F, {ptQAAxis, multQAAxis});

    // ---------------------- PID QA ----------------------

    registry.add("QA/PID/hTPCNSigma", "Distribution of TPC nSigma", kTH3F, {ptQAAxis, multQAAxis, nSigmaTPCQAaxis});
    registry.add("QA/PID/hTPCNSigmaK", "Distribution of TPC nSigma of K^{+} and K^{-}", kTH3F, {ptQAAxis, multQAAxis, nSigmaTPCQAaxis});
    registry.add("QA/PID/hTOFNSigma", "Distribution of TOF nSigma", kTH3F, {ptQAAxis, multQAAxis, nSigmaTOFQAaxis});
    registry.add("QA/PID/hTOFNSigmaK", "Distribution of TOF nSigma of K^{+} and K^{-}", kTH3F, {ptQAAxis, multQAAxis, nSigmaTOFQAaxis});
    registry.add("QA/PID/hTPCTOFnSigma", "", kTH3F, {ptQAAxis, nSigmaTPCQAaxis, nSigmaTOFQAaxis});

    registry.add("QA/PID/hTPCdEdxP", "dE/dx vs p of charged particles", kTH2F, {pQAaxis, dEdxQAaxis});
    registry.add("QA/PID/hTPCdEdxPK", "dE/dx vs p of K^{+} and K^{-}", kTH2F, {pQAaxis, dEdxQAaxis});
    registry.add("QA/PID/hTOFBetaP", "TOF #beta vs p of charged particles", kTH2F, {pQAaxis, betaQAaxis});
    registry.add("QA/PID/hTOFBetaPK", "TOF #beta vs p of K^{+} and K^{-}", kTH2F, {pQAaxis, betaQAaxis});

    // ------------------------- MC QA -------------------------
    if (static_cast<bool>(produce.produceMC)) {
      // Rec
      registry.add("QAMC/Rec/hSelection", "MC Rec True Event statistics", kTH1F, {{2, 0.0f, 2.0f}});
      auto hMCEventTruth = registry.get<TH1>(HIST("QAMC/Rec/hSelection"));
      hMCEventTruth->GetXaxis()->SetBinLabel(1, "Full MC Rec event statistics");
      hMCEventTruth->GetXaxis()->SetBinLabel(2, "MC Rec events passing event selection");
      hMCEventTruth->SetMinimum(0.1);

      // Gen
      registry.add("QAMC/Gen/hSelection", "MC Gen Event statistics", kTH1F, {{3, 0.0f, 3.0f}});
      auto hMCEventGen = registry.get<TH1>(HIST("QAMC/Gen/hSelection"));
      hMCEventGen->GetXaxis()->SetBinLabel(1, "Generated collisions");
      hMCEventGen->GetXaxis()->SetBinLabel(2, "Generated collisions with at least one reconstructed collision");
      hMCEventGen->GetXaxis()->SetBinLabel(3, "Generated collisions passing event selection");
      hMCEventGen->SetMinimum(0.1);

      // Factors
      registry.add("QAMC/Factors/hGenEvents", "Generated events", HistType::kTH2F, {nchQAAxis, {4, 0, 4}});
      auto hGenEvents = registry.get<TH2>(HIST("QAMC/Factors/hGenEvents"));
      hGenEvents->GetYaxis()->SetBinLabel(1, "All generated events");
      hGenEvents->GetYaxis()->SetBinLabel(2, "All reconstructed events");
      hGenEvents->GetYaxis()->SetBinLabel(3, "Generated events with at least one reconstructed event");
      hGenEvents->GetYaxis()->SetBinLabel(4, "Generated events passing event selection");

      registry.add("QAMC/Factors/hRecEvents", "Reconstructed events", HistType::kTH2F, {centQAAxis, {2, 0, 2}});
      auto hRecEvents = registry.get<TH2>(HIST("QAMC/Factors/hRecEvents"));
      hRecEvents->GetYaxis()->SetBinLabel(1, "All reconstructed events");
      hRecEvents->GetYaxis()->SetBinLabel(2, "Passing event selection");

      registry.add("QAMC/Factors/hGenALORESelEvents", "Centrality vs. Multiplicity of Generated Events with at least one reconstructed event passing event selection", kTH2F, {centQAAxis, nchQAAxis});
      registry.add("QAMC/Factors/hGenEventsCentNch", "Event centrality vs MC multiplicity", kTH2F, {centQAAxis, nchQAAxis});
      registry.add("QAMC/Factors/hNrecInGen", "Number of collisions in MC", kTH1F, {{10, -0.5, 9.5}});

      registry.add("QAMC/Factors/hGenPhi", "Generated #Phi", kTH3D, {nchQAAxis, centQAAxis, ptAxis});
      registry.add("QAMC/Factors/hGenALOREPhi", "Generated #Phi in collisions with at least one reconstructed collision", kTH3F, {nchQAAxis, centQAAxis, ptAxis});
      registry.add("QAMC/Factors/hRecPhi", "Reconstructed #Phi", kTH2F, {centQAAxis, ptAxis});

      // Resolution
      registry.add("QAMC/Resolution/h2ResolutionVz", "Resolution of collision V_{z}", kTH2F, {vzaxis, axisResolutionVz});
      auto hResVz = registry.get<TH2>(HIST("QAMC/Resolution/h2ResolutionVz"));
      hResVz->GetXaxis()->SetTitle("V_{z}^{rec} (cm)");
      hResVz->GetYaxis()->SetTitle("#DeltaV_{z} = V_{z}^{rec} - V_{z}^{gen} (cm)");

      registry.add("QAMC/Resolution/h2ResolutionPt", "Resolution of charged particles p_{T}", kTH2F, {ptQAAxis, axisResolutionPt});
      auto hResPt = registry.get<TH2>(HIST("QAMC/Resolution/h2ResolutionPt"));
      hResPt->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
      hResPt->GetYaxis()->SetTitle("#Deltap_{T} = p_{T}^{rec} - p_{T}^{gen} (GeV/c)");

      registry.add("QAMC/Resolution/h2ResolutionPtPhi", "p_{T} resolution vs p_{T}^{rec}", kTH2F, {ptQAAxis, axisResolutionPtPhi});
      auto hResPtPhi = registry.get<TH2>(HIST("QAMC/Resolution/h2ResolutionPtPhi"));
      hResPtPhi->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
      hResPtPhi->GetYaxis()->SetTitle("#Deltap_{T} = p_{T}^{rec} - p_{T}^{gen} (GeV/c)");

      registry.add("QAMC/Resolution/h2MassResolution", "Mass resolution vs p_{T}^{rec}", kTH2F, {ptQAAxis, axisResolutionMass});
      auto hResMass = registry.get<TH2>(HIST("QAMC/Resolution/h2MassResolution"));
      hResMass->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
      hResMass->GetYaxis()->SetTitle("#Deltam = m^{gen}_{KK} - m^{rec}_{KK} (GeV/c^{2})");
    }

    // ----------------------- Phi candidate QA -----------------------
    registry.add("QA/Phi/hRapidity", "Rapidity distribution of #Phi candidates", kTH3F, {ptQAAxis, multQAAxis, rapidityQAaxis});
    registry.add("QA/Phi/hEta", "Pseudorapidity distribution of #Phi candidates", kTH3F, {ptQAAxis, multQAAxis, etaQAaxis});
    registry.add("QA/Phi/hdPhi", "Azimuthal distribution (#Delta#phi) of #Phi candidates", kTH3F, {ptQAAxis, multQAAxis, dPhiQAaxis});
    registry.add("QA/Phi/hdPhideta", "Azimuthal distribution (#Delta#phi) of #Phi candidates vs #eta", kTH2F, {dEtaQAaxis, dPhiQAaxis});
    registry.add("QA/Phi/hdTheta", "Polar distribution (#Delta#theta) of #Phi candidates vs p_{T}", kTH3F, {ptQAAxis, multQAAxis, dThetaQAaxis});

    // Rotational background QA
    if (static_cast<bool>(produce.produceRotational)) {
      // Rotation around z axis
      registry.add("QA/RotationZ/hRapidity", "Rapidity distribution of #Phi candidates from rotational background", kTH3F, {ptQAAxis, multQAAxis, rapidityQAaxis});
      registry.add("QA/RotationZ/hEta", "Pseudorapidity distribution of #Phi candidates from rotational background", kTH3F, {ptQAAxis, multQAAxis, etaQAaxis});
      registry.add("QA/RotationZ/hdPhi", "Rotational background: Azimuthal distribution (#Delta#phi)", kTH3F, {ptQAAxis, multQAAxis, dPhiQAaxis});
      registry.add("QA/RotationZ/hdPhideta", "Rotational background: Azimuthal distribution (#Delta#phi) vs #eta", kTH2F, {dEtaQAaxis, dPhiQAaxis});
      registry.add("QA/RotationZ/hdTheta", "Rotational background: Polar distribution (#Delta#theta) vs p_{T}", kTH3F, {ptQAAxis, multQAAxis, dThetaQAaxis});
      // Momentum-axis rotation
      registry.add("QA/Rotation/hRapidity", "Rapidity distribution of #Phi candidates from rotational background", kTH3F, {ptQAAxis, multQAAxis, rapidityQAaxis});
      registry.add("QA/Rotation/hEta", "Pseudorapidity distribution of #Phi candidates from rotational background", kTH3F, {ptQAAxis, multQAAxis, etaQAaxis});
      registry.add("QA/Rotation/hdPhi", "Rotational background: Azimuthal distribution (#Delta#phi)", kTH3F, {ptQAAxis, multQAAxis, dPhiQAaxis});
      registry.add("QA/Rotation/hdPhideta", "Rotational background: Azimuthal distribution (#Delta#phi) vs #eta", kTH2F, {dEtaQAaxis, dPhiQAaxis});
      registry.add("QA/Rotation/hdTheta", "Rotational background: Polar distribution (#Delta#theta) vs p_{T}", kTH3F, {ptQAAxis, multQAAxis, dThetaQAaxis});
    }

    // Mixing QA
    if (mixingType != rsn::MixingType::none) {

      registry.add("QA/Mixing/h2mu1_mu2", "Event Mixing Multiplicity", kTH2F, {axisMultiplicityMixing, axisMultiplicityMixing});
      auto h2EMmu = registry.get<TH2>(HIST("QA/Mixing/h2mu1_mu2"));
      h2EMmu->GetXaxis()->SetTitle("1.Event multiplicity");
      h2EMmu->GetYaxis()->SetTitle("2.Event multiplicity");

      registry.add("QA/Mixing/h2ce1_ce2", "Event Mixing Centrality", kTH2F, {axisCentralityMixing, axisCentralityMixing});
      auto h2EMce = registry.get<TH2>(HIST("QA/Mixing/h2ce1_ce2"));
      h2EMce->GetXaxis()->SetTitle("1.Event centrality");
      h2EMce->GetYaxis()->SetTitle("2.Event centrality");

      registry.add("QA/Mixing/h2vz1_vz2", "Event Mixing Vertex z", kTH2F, {axisVertexMixing, axisVertexMixing});
      auto hEMTvz = registry.get<TH2>(HIST("QA/Mixing/h2vz1_vz2"));
      hEMTvz->GetXaxis()->SetTitle("1.Event V_{z}");
      hEMTvz->GetYaxis()->SetTitle("2.Event V_{z}");

      registry.add("QA/Mixing/hdPhideta", "Mixing background: Azimuthal distribution (#Delta#phi) vs #eta", kTH2F, {dEtaQAaxis, dPhiQAaxis});
    }
  }
  template <typename T>
  bool selectedEvent(const T& collision)
  {
    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 0.5); // all events
    }

    if (static_cast<bool>(eventCuts.isTriggerTVX) && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 1.5); // events passing trigger TVX cut
    }

    if (static_cast<bool>(eventCuts.noTimeFrameBorder) && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 2.5); // events passing no time frame border cut
    }

    if (static_cast<bool>(eventCuts.noITSROFrameBorder) && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 3.5); // events passing no ITS RO frame border cut
    }

    if (static_cast<bool>(eventCuts.sel8) && !collision.sel8()) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 4.5); // events passing sel8 cut (contains all the previous cuts)
    }

    if (static_cast<bool>(eventCuts.isVertexITSTPC) && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 5.5); // events passing IsVertexITSTPC cut
    }

    if (static_cast<bool>(eventCuts.noSameBunchPileup) && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 6.5); // events passing no same bunch pileup cut
    }

    if (static_cast<bool>(eventCuts.isGoodZvtxFT0vsPV) && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 7.5); // events passing IsGoodZvtxFT0vsPV cut
    }

    if (std::abs(collision.posZ()) > static_cast<float>(eventCuts.vzCut)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 8.5); // events passing V_{z} cut
    }

    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 9.5); // INEL
    }

    if (static_cast<bool>(eventCuts.inelGt0) && !collision.isInelGt0()) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Event/hSelection"), 10.5); // events passing INEL>0 cut
      registry.fill(HIST("QA/Event/hVtxZ"), collision.posZ());
      registry.fill(HIST("QA/Event/hMult"), getMultiplicity(collision));
      registry.fill(HIST("QA/Event/hCent"), getCentrality(collision));
    }

    return true;
  }
  template <typename T>
  float tpcNsigma(const T& track)
  {
    float tpcNsigma = 0.0f;
    int particleType = (track.sign() > 0) ? static_cast<int>(daughterPos) : static_cast<int>(daughterNeg);

    if (particleType == pion) {
      tpcNsigma = track.tpcNSigmaPi();
    } else if (particleType == kaon) {
      tpcNsigma = track.tpcNSigmaKa();
    } else if (particleType == proton) {
      tpcNsigma = track.tpcNSigmaPr();
    }
    return tpcNsigma;
  }
  template <typename T>
  float tofNsigma(const T& track)
  {
    float tofNsigma = 0.0f;
    int particleType = (track.sign() > 0) ? static_cast<int>(daughterPos) : static_cast<int>(daughterNeg);

    if (particleType == pion) {
      tofNsigma = track.tofNSigmaPi();
    } else if (particleType == kaon) {
      tofNsigma = track.tofNSigmaKa();
    } else if (particleType == proton) {
      tofNsigma = track.tofNSigmaPr();
    }
    return tofNsigma;
  }
  template <typename T>
  bool selectedTrack(const T& track, bool isPositive)
  {
    if (dataQA) {
      registry.fill(HIST("QA/Track/hSelection"), 0.5);
    }

    if (track.pt() < static_cast<float>(trackCuts.pt)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Track/hSelection"), 1.5);
    }

    if (std::abs(track.eta()) >= static_cast<float>(trackCuts.etatrack)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Track/hSelection"), 2.5);
    }

    if (std::abs(track.dcaXY()) >= static_cast<float>(trackCuts.dcaXY) ||
        std::abs(track.dcaZ()) >= static_cast<float>(trackCuts.dcaZ)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Track/hSelection"), 3.5);
    }

    // PID selection: TPC-only for pt < threshold value, TPC+TOF for pt >= threshold value and have TOF, else TPC-only
    float nSigmaCut = isPositive ? static_cast<float>(trackCuts.tpcnSigmaPos) : static_cast<float>(trackCuts.tpcnSigmaNeg);
    if (track.pt() < static_cast<float>(trackCuts.ptTOFThreshold) || !track.hasTOF() || static_cast<bool>(trackCuts.tpcPidOnly)) {
      if (std::abs(tpcNsigma(track)) >= nSigmaCut) {
        return false;
      }
    } else {
      if (std::sqrt(tpcNsigma(track) * tpcNsigma(track) + tofNsigma(track) * tofNsigma(track)) >= static_cast<float>(trackCuts.combinedNSigma)) {
        return false;
      }
    }
    if (dataQA) {
      registry.fill(HIST("QA/Track/hSelection"), 4.5);
    }

    if (track.tpcNClsFound() < static_cast<int>(trackCuts.tpcNClsFound)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Track/hSelection"), 5.5);
    }

    if (track.tpcNClsCrossedRows() < static_cast<int>(trackCuts.tpcNClsCrossedRows)) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Track/hSelection"), 6.5);
    }

    if (static_cast<bool>(trackCuts.globalTrack)) {
      if (!track.isGlobalTrack()) {
        return false;
      }
    } else if (static_cast<bool>(trackCuts.primaryTrack)) {
      if (!track.isPrimaryTrack()) {
        return false;
      }
    }
    if (dataQA) {
      registry.fill(HIST("QA/Track/hSelection"), 7.5);
    }

    if (static_cast<bool>(trackCuts.pvContributor) && !track.isPVContributor()) {
      return false;
    }
    if (dataQA) {
      registry.fill(HIST("QA/Track/hSelection"), 8.5);
    }

    return true;
  }
  template <typename T>
  ROOT::Math::PxPyPzMVector calculateMother(const T& track1, const T& track2)
  {
    d1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPos);
    d2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massNeg);
    return d1 + d2;
  }
  bool selectedMother(const ROOT::Math::PxPyPzMVector& motherCandidate)
  {
    return std::abs(motherCandidate.Rapidity()) <= static_cast<float>(trackCuts.rapidity);
  }
  template <typename T>
  float getMultiplicity(const T& collision)
  {
    float multiplicity = collision.multFT0M();
    return multiplicity;
  }
  template <typename T>
  float getCentrality(const T& collision)
  {
    float centrality = collision.centFT0M();
    return centrality;
  }
  double* fillPointPair(double im, double pt, double mu, double ce, double ns1, double ns2, double eta, double y, double vz, double mum, double cem, double vzm)
  {
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::im)] = im;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::pt)] = pt;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::mu)] = mu;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ce)] = ce;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns1)] = ns1;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns2)] = ns2;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::eta)] = eta;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::y)] = y;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::vz)] = vz;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::mum)] = mum;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::cem)] = cem;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::vzm)] = vzm;

    return pointPair;
  }

  void processQA(EventCandidate const& collision, TrackCandidates const& tracks)
  {
    dataQA = true;
    bool selected = selectedEvent(collision);
    dataQA = false;

    if (!selected) {
      return;
    }

    double centrality = getCentrality(collision);

    int nch = 0;
    for (const auto& track : tracks) {

      registry.fill(HIST("QA/Track/hEta"), track.pt(), centrality, track.eta());
      registry.fill(HIST("QA/Track/hPt"), track.pt(), centrality);
      registry.fill(HIST("QA/Track/hDCAxy"), track.pt(), centrality, track.dcaXY());
      registry.fill(HIST("QA/Track/hDCAz"), track.pt(), centrality, track.dcaZ());
      registry.fill(HIST("QA/Track/hTPCNClsFound"), track.pt(), centrality, track.tpcNClsFound());
      registry.fill(HIST("QA/Track/hTPCNClsCrossedRows"), track.pt(), centrality, track.tpcNClsCrossedRows());
      registry.fill(HIST("QA/Track/hRapidity"), track.pt(), centrality, track.sign() > 0 ? track.rapidity(massPos) : track.rapidity(massNeg));

      registry.fill(HIST("QA/PID/hTPCNSigma"), track.pt(), centrality, tpcNsigma(track));
      registry.fill(HIST("QA/PID/hTPCdEdxP"), track.p(), track.tpcSignal());

      if (track.hasTOF()) {
        registry.fill(HIST("QA/PID/hTOFNSigma"), track.pt(), centrality, tofNsigma(track));
        registry.fill(HIST("QA/PID/hTOFBetaP"), track.p(), track.beta());
      }

      if (track.isPrimaryTrack() && std::abs(track.eta()) < static_cast<float>(trackCuts.etatrack)) {
        nch++;
      }

      dataQA = true;
      bool selectedTrackCandidate = selectedTrack(track, track.sign() > 0);
      dataQA = false;

      if (!selectedTrackCandidate) {
        continue;
      }

      registry.fill(HIST("QA/Kaon/hEta"), track.pt(), centrality, track.eta());
      registry.fill(HIST("QA/Kaon/hPt"), track.pt(), centrality);
      registry.fill(HIST("QA/Kaon/hDCAxy"), track.pt(), centrality, track.dcaXY());
      registry.fill(HIST("QA/Kaon/hDCAz"), track.pt(), centrality, track.dcaZ());
      registry.fill(HIST("QA/Kaon/hTPCNClsFound"), track.pt(), centrality, track.tpcNClsFound());
      registry.fill(HIST("QA/Kaon/hTPCNClsCrossedRows"), track.pt(), centrality, track.tpcNClsCrossedRows());
      registry.fill(HIST("QA/Kaon/hRapidity"), track.pt(), centrality, track.sign() > 0 ? track.rapidity(massPos) : track.rapidity(massNeg));

      registry.fill(HIST("QA/PID/hTPCNSigmaK"), track.pt(), centrality, tpcNsigma(track));
      registry.fill(HIST("QA/PID/hTPCdEdxPK"), track.p(), track.tpcSignal());

      if (track.hasTOF()) {
        registry.fill(HIST("QA/PID/hTOFNSigmaK"), track.pt(), centrality, tofNsigma(track));
        registry.fill(HIST("QA/PID/hTOFBetaPK"), track.p(), track.beta());
        registry.fill(HIST("QA/PID/hTPCTOFnSigma"), track.pt(), tpcNsigma(track), tofNsigma(track));
      }
    }
    registry.fill(HIST("QA/Event/hCentNch"), getCentrality(collision), nch);
  }
  PROCESS_SWITCH(PhianalysisTHnSparse, processQA, "Process Event for Data", true);

  void processData(EventCandidate const& collision, TrackCandidates const& /*tracks*/)
  {
    auto posDaughters = positive->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDaughters = negative->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (!selectedEvent(collision)) {
      return;
    }

    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughters, negDaughters))) {

      if (!selectedTrack(track1, true)) {
        continue;
      }
      if (!selectedTrack(track2, false)) {
        continue;
      }

      mother = calculateMother(track1, track2);
      if (!selectedMother(mother)) {
        continue;
      }

      registry.fill(HIST("QA/Phi/hRapidity"), mother.Pt(), getCentrality(collision), mother.Rapidity());
      registry.fill(HIST("QA/Phi/hEta"), mother.Pt(), getCentrality(collision), mother.Eta());
      registry.fill(HIST("QA/Phi/hdPhi"), mother.Pt(), getCentrality(collision), track1.phi() - track2.phi());
      registry.fill(HIST("QA/Phi/hdPhideta"), track1.eta() - track2.eta(), track1.phi() - track2.phi());
      registry.fill(HIST("QA/Phi/hdTheta"), mother.Pt(), getCentrality(collision), d1.Theta() - d2.Theta());

      pointPair = fillPointPair(mother.M(),
                                mother.Pt(),
                                getMultiplicity(collision),
                                getCentrality(collision),
                                tpcNsigma(track1),
                                tpcNsigma(track2),
                                mother.Eta(),
                                mother.Rapidity(),
                                collision.posZ(),
                                0,
                                0,
                                0);
      rsnOutput->fillUnlikepm(pointPair);

      if (static_cast<bool>(produce.produceRotational)) {

        // Rotational background with rotation of track2 around track1
        for (int i = 1; i <= static_cast<int>(nRotations); i++) {
          float starting = static_cast<float>(startingAngle) * TMath::DegToRad();
          float angle = starting + i * ((o2::constants::math::TwoPI - 2 * starting) / (static_cast<int>(nRotations) + 1));
          float px2new = track2.px() * std::cos(angle) - track2.py() * std::sin(angle);
          float py2new = track2.px() * std::sin(angle) + track2.py() * std::cos(angle);
          ROOT::Math::PxPyPzMVector d2rot(px2new, py2new, track2.pz(), massNeg);
          auto motherRotZ = d1 + d2rot;

          registry.fill(HIST("QA/RotationZ/hRapidity"), motherRotZ.Pt(), getCentrality(collision), motherRotZ.Rapidity());
          registry.fill(HIST("QA/RotationZ/hEta"), motherRotZ.Pt(), getCentrality(collision), motherRotZ.Eta());
          registry.fill(HIST("QA/RotationZ/hdPhi"), motherRotZ.Pt(), getCentrality(collision), d1.Phi() - d2rot.Phi());
          registry.fill(HIST("QA/RotationZ/hdPhideta"), d1.Eta() - d2rot.Eta(), d1.Phi() - d2rot.Phi());
          registry.fill(HIST("QA/RotationZ/hdTheta"), motherRotZ.Pt(), getCentrality(collision), d1.Theta() - d2rot.Theta());

          pointPair = fillPointPair(motherRotZ.M(),
                                    motherRotZ.Pt(),
                                    getMultiplicity(collision),
                                    getCentrality(collision),
                                    tpcNsigma(track1),
                                    tpcNsigma(track2),
                                    motherRotZ.Eta(),
                                    motherRotZ.Rapidity(),
                                    collision.posZ(),
                                    0,
                                    0,
                                    0);

          rsnOutput->fillRotationZ(pointPair);
        }

        // Rotational background with rotation of 90 degrees around mother momentum axis
        ROOT::Math::AxisAngle rotationAxis(mother.Vect(), constants::math::PIHalf);
        ROOT::Math::Rotation3D rotationMatrix(rotationAxis);

        const auto rotD1 = rotationMatrix * d1;
        const auto rotD2 = rotationMatrix * d2;

        if (negDaughters.size() > 1) {
          for (const auto& track3 : negDaughters) {
            if (track3.globalIndex() == track1.globalIndex() || track3.globalIndex() == track2.globalIndex()) {
              continue;
            }
            if (!selectedTrack(track3, false)) {
              continue;
            }
            ROOT::Math::PxPyPzMVector d3(track3.px(), track3.py(), track3.pz(), massNeg);
            auto motherRot = rotD1 + d3;

            registry.fill(HIST("QA/Rotation/hRapidity"), motherRot.Pt(), getCentrality(collision), motherRot.Rapidity());
            registry.fill(HIST("QA/Rotation/hEta"), motherRot.Pt(), getCentrality(collision), motherRot.Eta());
            registry.fill(HIST("QA/Rotation/hdPhi"), motherRot.Pt(), getCentrality(collision), rotD1.Phi() - d3.Phi());
            registry.fill(HIST("QA/Rotation/hdPhideta"), rotD1.Eta() - d3.Eta(), rotD1.Phi() - d3.Phi());
            registry.fill(HIST("QA/Rotation/hdTheta"), motherRot.Pt(), getCentrality(collision), rotD1.Theta() - d3.Theta());

            pointPair = fillPointPair(motherRot.M(),
                                      motherRot.Pt(),
                                      getMultiplicity(collision),
                                      getCentrality(collision),
                                      tpcNsigma(track1),
                                      tpcNsigma(track2),
                                      motherRot.Eta(),
                                      motherRot.Rapidity(),
                                      collision.posZ(),
                                      0,
                                      0,
                                      0);

            rsnOutput->fillRotation(pointPair);

            // Likesign rotation for negative daughters
            motherRot = rotD2 + d3;
            pointPair = fillPointPair(motherRot.M(),
                                      motherRot.Pt(),
                                      getMultiplicity(collision),
                                      getCentrality(collision),
                                      tpcNsigma(track1),
                                      tpcNsigma(track2),
                                      motherRot.Eta(),
                                      motherRot.Rapidity(),
                                      collision.posZ(),
                                      0,
                                      0,
                                      0);

            rsnOutput->fillRotationLike(pointPair);
          }
        }

        if (posDaughters.size() > 1) {
          for (const auto& track3 : posDaughters) {

            if (track3.globalIndex() == track1.globalIndex() || track3.globalIndex() == track2.globalIndex()) {
              continue;
            }

            if (!selectedTrack(track3, true)) {
              continue;
            }

            ROOT::Math::PxPyPzMVector d3(track3.px(), track3.py(), track3.pz(), massPos);

            auto motherRot = rotD2 + d3;

            registry.fill(HIST("QA/Rotation/hRapidity"), motherRot.Pt(), getCentrality(collision), motherRot.Rapidity());
            registry.fill(HIST("QA/Rotation/hEta"), motherRot.Pt(), getCentrality(collision), motherRot.Eta());
            registry.fill(HIST("QA/Rotation/hdPhi"), motherRot.Pt(), getCentrality(collision), rotD2.Phi() - d3.Phi());
            registry.fill(HIST("QA/Rotation/hdPhideta"), rotD2.Eta() - d3.Eta(), rotD2.Phi() - d3.Phi());
            registry.fill(HIST("QA/Rotation/hdTheta"), motherRot.Pt(), getCentrality(collision), rotD2.Theta() - d3.Theta());

            pointPair = fillPointPair(motherRot.M(),
                                      motherRot.Pt(),
                                      getMultiplicity(collision),
                                      getCentrality(collision),
                                      tpcNsigma(track1),
                                      tpcNsigma(track2),
                                      motherRot.Eta(),
                                      motherRot.Rapidity(),
                                      collision.posZ(),
                                      0,
                                      0,
                                      0);

            rsnOutput->fillRotation(pointPair);

            // Likesign rotation for positive daughters
            motherRot = rotD1 + d3;
            pointPair = fillPointPair(motherRot.M(),
                                      motherRot.Pt(),
                                      getMultiplicity(collision),
                                      getCentrality(collision),
                                      tpcNsigma(track1),
                                      tpcNsigma(track2),
                                      motherRot.Eta(),
                                      motherRot.Rapidity(),
                                      collision.posZ(),
                                      0,
                                      0,
                                      0);

            rsnOutput->fillRotationLike(pointPair);
          }
        }
      }
    }
    if (static_cast<bool>(produce.produceLikesign)) {

      for (const auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posDaughters, posDaughters))) {
        if (!selectedTrack(track1, true)) {
          continue;
        }
        if (!selectedTrack(track2, true)) {
          continue;
        }

        mother = calculateMother(track1, track2);
        if (!selectedMother(mother)) {
          continue;
        }

        pointPair = fillPointPair(mother.M(),
                                  mother.Pt(),
                                  getMultiplicity(collision),
                                  getCentrality(collision),
                                  tpcNsigma(track1),
                                  tpcNsigma(track2),
                                  mother.Eta(),
                                  mother.Rapidity(),
                                  collision.posZ(),
                                  0,
                                  0,
                                  0);

        rsnOutput->fillLikepp(pointPair);
      }

      for (const auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(negDaughters, negDaughters))) {
        if (!selectedTrack(track1, false)) {
          continue;
        }
        if (!selectedTrack(track2, false)) {
          continue;
        }

        mother = calculateMother(track1, track2);
        if (!selectedMother(mother)) {
          continue;
        }

        pointPair = fillPointPair(mother.M(),
                                  mother.Pt(),
                                  getMultiplicity(collision),
                                  getCentrality(collision),
                                  tpcNsigma(track1),
                                  tpcNsigma(track2),
                                  mother.Eta(),
                                  mother.Rapidity(),
                                  collision.posZ(),
                                  0,
                                  0,
                                  0);

        rsnOutput->fillLikemm(pointPair);
      }
    }
  }
  PROCESS_SWITCH(PhianalysisTHnSparse, processData, "Process Event for Data", true);

  void processTrue(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    if (!static_cast<bool>(produce.produceMC)) {
      return;
    }

    registry.fill(HIST("QAMC/Rec/hSelection"), 0.5);
    registry.fill(HIST("QAMC/Factors/hRecEvents"), getCentrality(collision), 0.5);

    if (!selectedEvent(collision)) {
      return;
    }

    registry.fill(HIST("QAMC/Rec/hSelection"), 1.5);
    registry.fill(HIST("QAMC/Factors/hRecEvents"), getCentrality(collision), 1.5);

    auto posDaughtersMC = positiveMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDaughtersMC = negativeMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (!collision.has_mcCollision()) {
      return;
    }

    auto mcCollision = collision.mcCollision();
    registry.fill(HIST("QAMC/Resolution/h2ResolutionVz"), collision.posZ(), (collision.posZ() - mcCollision.posZ()));

    for (const auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto mctrack = track.mcParticle();
        registry.fill(HIST("QAMC/Resolution/h2ResolutionPt"), track.pt(), (track.pt() - mctrack.pt()));
      }
    }

    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersMC, negDaughtersMC))) {

      if (!track1.has_mcParticle()) {
        continue;
      }
      if (!track2.has_mcParticle()) {
        continue;
      }

      if (!selectedTrack(track1, true)) {
        continue;
      }
      if (!selectedTrack(track2, false)) {
        continue;
      }

      const auto mctrack1 = track1.mcParticle();
      const auto mctrack2 = track2.mcParticle();
      int track1PDG = std::abs(mctrack1.pdgCode());
      int track2PDG = std::abs(mctrack2.pdgCode());

      if (track1PDG != daughterPosPDG || track2PDG != daughterNegPDG) {
        continue;
      }

      int n = 0;
      for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
        for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {

          if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
            continue;
          }

          if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
            continue;
          }

          if (std::abs(mothertrack1.y()) > static_cast<float>(trackCuts.rapidity)) {
            continue;
          }

          if (std::abs(mothertrack1.pdgCode()) != motherPDG) {
            continue;
          }

          mother = calculateMother(track1, track2);
          motherGen = calculateMother(mctrack1, mctrack2);
          if (!selectedMother(mother)) {
            continue;
          }

          if (n > 0) {
            continue;
          }

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    getMultiplicity(collision),
                                    getCentrality(collision),
                                    tpcNsigma(track1),
                                    tpcNsigma(track2),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    collision.posZ(),
                                    0,
                                    0,
                                    0);

          registry.fill(HIST("QAMC/Resolution/h2ResolutionPtPhi"), mother.Pt(), (mother.Pt() - mothertrack1.pt()));
          registry.fill(HIST("QAMC/Resolution/h2MassResolution"), mother.Pt(), (motherGen.M() - mother.M()));

          rsnOutput->fillUnlikeTrueRec(pointPair);

          pointPair = fillPointPair(motherGen.M(),
                                    motherGen.Pt(),
                                    getMultiplicity(collision),
                                    getCentrality(collision),
                                    tpcNsigma(track1),
                                    tpcNsigma(track2),
                                    motherGen.Eta(),
                                    motherGen.Rapidity(),
                                    collision.posZ(),
                                    0,
                                    0,
                                    0);

          rsnOutput->fillUnlikeTrueGen(pointPair);

          registry.fill(HIST("QAMC/Factors/hRecPhi"), getCentrality(collision), motherGen.Pt());

          n++;
        }
      }
    }
  }
  PROCESS_SWITCH(PhianalysisTHnSparse, processTrue, "Process Event for MC reconstruction.", false);

  void processGen(McCollisionMults::iterator const& mcCollision, soa::SmallGroups<EventCandidatesMCGen> const& collisions, LabeledTracks const& /*particles*/, aod::McParticles const& mcParticles)
  {
    if (!static_cast<bool>(produce.produceMC)) {
      return;
    }

    registry.fill(HIST("QAMC/Factors/hNrecInGen"), collisions.size());
    registry.fill(HIST("QAMC/Gen/hSelection"), 0.5);
    registry.fill(HIST("QAMC/Factors/hGenEvents"), mcCollision.multMCNParticlesEta05(), 0.5);

    if (collisions.size() > 0) {
      registry.fill(HIST("QAMC/Gen/hSelection"), 1.5);
      registry.fill(HIST("QAMC/Factors/hGenEvents"), mcCollision.multMCNParticlesEta05(), 2.5);
    }

    int nContributors = -1;
    bool hasSelectedCollision = false;
    float centrality = 100.5f;
    float multiplicity = 0.f;

    for (const auto& collision : collisions) {
      registry.fill(HIST("QAMC/Factors/hGenEvents"), mcCollision.multMCNParticlesEta05(), 1.5);
      if (!selectedEvent(collision)) {
        continue;
      }

      if (collision.numContrib() > nContributors) {
        nContributors = collision.numContrib();
        centrality = getCentrality(collision);
        multiplicity = getMultiplicity(collision);
        hasSelectedCollision = true;
      }
    }

    registry.fill(HIST("QAMC/Factors/hGenEventsCentNch"), centrality, mcCollision.multMCNParticlesEta05());

    // All generated Phi mesons
    for (const auto& particle : mcParticles) {

      if (std::abs(particle.y()) > static_cast<float>(trackCuts.rapidity)) {
        continue;
      }

      if (particle.pdgCode() == motherPDG) {

        auto daughters = particle.daughters_as<aod::McParticles>();
        if (daughters.size() != dauSize) {
          continue;
        }

        auto daup = false;
        auto daun = false;

        for (const auto& dau : daughters) {
          if (dau.pdgCode() == daughterPosPDG) {
            daup = true;
            d1 = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), massPos);
          } else if (dau.pdgCode() == -daughterNegPDG) {
            daun = true;
            d2 = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), massNeg);
          }
        }
        if (!daup || !daun) {
          continue;
        }

        mother = d1 + d2;

        registry.fill(HIST("QAMC/Factors/hGenPhi"), mcCollision.multMCNParticlesEta05(), centrality, particle.pt());
      }
    }

    if (!hasSelectedCollision) {
      return;
    }

    registry.fill(HIST("QAMC/Gen/hSelection"), 2.5);
    registry.fill(HIST("QAMC/Factors/hGenEvents"), mcCollision.multMCNParticlesEta05(), 3.5);
    registry.fill(HIST("QAMC/Factors/hGenALORESelEvents"), centrality, mcCollision.multMCNParticlesEta05());

    // Generated Phi mesons in selected collisions
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) > static_cast<float>(trackCuts.rapidity)) {
        continue;
      }

      if (mcParticle.pdgCode() == motherPDG) {
        auto daughters = mcParticle.daughters_as<aod::McParticles>();
        if (daughters.size() != dauSize) {
          continue;
        }

        auto daup = false;
        auto daun = false;

        for (const auto& dau : daughters) {
          if (dau.pdgCode() == daughterPosPDG) {
            daup = true;
            d1 = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), massPos);
          } else if (dau.pdgCode() == -daughterNegPDG) {
            daun = true;
            d2 = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), massNeg);
          }
        }

        if (!daup || !daun) {
          continue;
        }

        mother = d1 + d2;

        pointPair = fillPointPair(mother.M(),
                                  mother.Pt(),
                                  multiplicity,
                                  centrality,
                                  0,
                                  0,
                                  mother.Eta(),
                                  mother.Rapidity(),
                                  mcCollision.posZ(),
                                  0,
                                  0,
                                  0);

        rsnOutput->fillUnlikeGen(pointPair);

        registry.fill(HIST("QAMC/Factors/hGenALOREPhi"), mcCollision.multMCNParticlesEta05(), centrality, mother.Pt());
      }
    }
  }

  PROCESS_SWITCH(PhianalysisTHnSparse, processGen, "Process MC Generated.", false);

  void processMixed(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    if (mixingType == rsn::MixingType::none) {
      return;
    }

    auto tracksTuple = std::make_tuple(tracks);

    BinningTypeVzCe binningVzCe{{axisVertexMixing, axisCentralityMixing}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVzCe> pairVzCe{binningVzCe, static_cast<int>(nMixedEvents), -1, collisions, tracksTuple, &cache};

    BinningTypeVzMu binningVzMu{{axisVertexMixing, axisMultiplicityMixing}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVzMu> pairVzMu{binningVzMu, static_cast<int>(nMixedEvents), -1, collisions, tracksTuple, &cache};

    if (mixingType == rsn::MixingType::ce) {
      for (const auto& [c1, tracks1, c2, tracks2] : pairVzCe) {
        if (!selectedEvent(c1) || !selectedEvent(c2)) {
          continue;
        }

        auto posDaughtersc1 = positive->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto posDaughtersc2 = positive->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);
        auto negDaughtersc1 = negative->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto negDaughtersc2 = negative->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);

        registry.fill(HIST("QA/Mixing/h2mu1_mu2"), getMultiplicity(c1), getMultiplicity(c2));
        registry.fill(HIST("QA/Mixing/h2ce1_ce2"), getCentrality(c1), getCentrality(c2));
        registry.fill(HIST("QA/Mixing/h2vz1_vz2"), c1.posZ(), c2.posZ());

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersc1, negDaughtersc2))) {

          if (!selectedTrack(track1, true)) {
            continue;
          }
          if (!selectedTrack(track2, false)) {
            continue;
          }

          mother = calculateMother(track1, track2);
          if (!selectedMother(mother)) {
            continue;
          }

          registry.fill(HIST("QA/Mixing/hdPhideta"), track1.eta() - track2.eta(), track1.phi() - track2.phi());

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    tpcNsigma(track1),
                                    tpcNsigma(track2),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    c1.posZ(),
                                    getMultiplicity(c2),
                                    getCentrality(c2),
                                    c2.posZ());

          rsnOutput->fillMixingpm(pointPair);
        }

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersc2, negDaughtersc1))) {

          if (!selectedTrack(track1, true)) {
            continue;
          }
          if (!selectedTrack(track2, false)) {
            continue;
          }

          mother = calculateMother(track1, track2);
          if (!selectedMother(mother)) {
            continue;
          }

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    tpcNsigma(track1),
                                    tpcNsigma(track2),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    c1.posZ(),
                                    getMultiplicity(c2),
                                    getCentrality(c2),
                                    c2.posZ());

          rsnOutput->fillMixingmp(pointPair);
        }
      }
    }
    if (mixingType == rsn::MixingType::mu) {
      for (const auto& [c1, tracks1, c2, tracks2] : pairVzMu) {
        if (!selectedEvent(c1) || !selectedEvent(c2)) {
          continue;
        }

        auto posDaughtersc1 = positive->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto posDaughtersc2 = positive->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);
        auto negDaughtersc1 = negative->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto negDaughtersc2 = negative->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);

        registry.fill(HIST("QA/Mixing/h2mu1_mu2"), getMultiplicity(c1), getMultiplicity(c2));
        registry.fill(HIST("QA/Mixing/h2ce1_ce2"), getCentrality(c1), getCentrality(c2));
        registry.fill(HIST("QA/Mixing/h2vz1_vz2"), c1.posZ(), c2.posZ());

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersc1, negDaughtersc2))) {

          if (!selectedTrack(track1, true)) {
            continue;
          }

          if (!selectedTrack(track2, false)) {
            continue;
          }

          mother = calculateMother(track1, track2);
          if (!selectedMother(mother)) {
            continue;
          }

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    tpcNsigma(track1),
                                    tpcNsigma(track2),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    c1.posZ(),
                                    getMultiplicity(c2),
                                    getCentrality(c2),
                                    c2.posZ());

          rsnOutput->fillMixingpm(pointPair);
        }

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersc2, negDaughtersc1))) {

          if (!selectedTrack(track1, true)) {
            continue;
          }
          if (!selectedTrack(track2, false)) {
            continue;
          }

          mother = calculateMother(track1, track2);
          if (!selectedMother(mother)) {
            continue;
          }

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    tpcNsigma(track1),
                                    tpcNsigma(track2),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    c1.posZ(),
                                    getMultiplicity(c2),
                                    getCentrality(c2),
                                    c2.posZ());

          rsnOutput->fillMixingmp(pointPair);
        }
      }
    }
  }
  PROCESS_SWITCH(PhianalysisTHnSparse, processMixed, "Process Mixing Event.", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PhianalysisTHnSparse>(cfgc)};
}
