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

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

#include <Math/Vector4D.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PhianalysisTHnSparse {

  SliceCache cache;

  struct : ConfigurableGroup {
    Configurable<bool> produceQA{"produceQA", false, "Produce qa histograms."};
    Configurable<bool> produceStats{"produceStats", false, "Produce statistics histograms."};
    Configurable<bool> produceTrue{"produceTrue", false, "Produce True and Gen histograms."};
    Configurable<bool> produceLikesign{"produceLikesign", false, "Produce Like sign histograms."};
    Configurable<std::string> eventMixing{"eventMixing", "none", "Produce Event Mixing histograms of type."};
    Configurable<bool> produceRotational{"produceRotational", false, "Produce Rotational histograms."};
  } produce;

  Configurable<int> daughterPos{"daughterPos", 3, "Particle type of the positive dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> daughterNeg{"daughterNeg", 3, "Particle type of the negative dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> motherPDG{"motherPDG", 333, "PDG code of mother particle."};
  Configurable<int> daughterPosPDG{"daughterPosPDG", 321, "PDG code of positive dauther particle."};
  Configurable<int> daughterNegPDG{"daughterNegPDG", 321, "PDG code of negative dauther particle."};

  struct : ConfigurableGroup {
    Configurable<float> tpcnSigmaPos{"tpcnSigmaPos", 3.0f, "TPC NSigma cut of the positive particle."};
    Configurable<float> tpcnSigmaNeg{"tpcnSigmaNeg", 3.0f, "TPC NSigma cut of the negative particle."};
    Configurable<bool> tpcPidOnly{"tpcPidOnly", false, "Use TPC only for PID."};
    Configurable<float> combinedNSigma{"combinedNSigma", 3.0f, "Combined NSigma cut for combined TPC and TOF NSigma cut."};
    Configurable<float> ptTOFThreshold{"ptTOFThreshold", 0.5f, "Threshold for applying TOF."};
    Configurable<float> rapidity{"rapidity", 0.5f, "Rapidity cut (maximum)."};
    Configurable<float> etatrack{"etatrack", 0.8f, "Eta cut for track."};
    Configurable<float> pt{"pt", 0.15f, "Cut: Minimal value of tracks pt."};
    Configurable<float> dcaXY{"dcaXY", 1.0f, "Cut: Maximal value of tracks DCA XY."};
    Configurable<float> dcaZ{"dcaZ", 1.0f, "Cut: Maximal value of tracks DCA Z."};
    Configurable<bool> globalTrack{"globalTrack", false, "Use global track selection."};
    Configurable<bool> inelGrater0{"inelGrater0", true, "Select events with INEL>0."};
    Configurable<int> tpcNClsFound{"tpcNClsFound", 70, "Cut: Minimal value of found TPC clasters"};
    Configurable<float> vzCut{"vzCut", 10.0f, "Cut: Maximal value of Z vertex position."};
  } cut;

  struct : ConfigurableGroup {
    Configurable<int> verboselevel{"verboselevel", 0, "Verbose level"};
    Configurable<int> refresh{"refresh", 0, "Freqency of print event information."};
    Configurable<int> refreshIndex{"refreshIndex", 0, "Freqency of print event information index."};
  } verbose;

  Configurable<std::vector<std::string>> sparseAxes{"sparseAxes", std::vector<std::string>{o2::analysis::rsn::pair_axis::names}, "Axes."};
  Configurable<std::vector<std::string>> sysAxes{"sysAxes", std::vector<std::string>{o2::analysis::rsn::systematic_axis::names}, "Axes."};

  ConfigurableAxis invaxis{"invaxis", {130, 0.97, 1.1}, "Invariant mass axis binning."};
  ConfigurableAxis ptaxis{"ptaxis", {20, 0., 20.}, "Pt axis binning."};
  ConfigurableAxis vzaxis{"vzaxis", {40, -20., 20.}, "Z vertex position axis binning."};
  ConfigurableAxis multiplicityaxis{"multiplicityaxis", {50, 0., 5000.}, "Multiplicity axis binning."};
  ConfigurableAxis centralityaxis{"centralityaxis", {20, 0., 100.}, "Centrality axis binning."};
  ConfigurableAxis etaaxis{"etaaxis", {16., -1.0 * static_cast<float>(cut.etatrack), static_cast<float>(cut.etatrack)}, "Pseudorapidity axis binning."};
  ConfigurableAxis rapidityaxis{"rapidityaxis", {10., -1.0 * static_cast<float>(cut.rapidity), static_cast<float>(cut.rapidity)}, "Rapidity axis binning."};
  ConfigurableAxis nsigmaaxisPos{"nsigmaaxisPos", {1, -static_cast<float>(cut.tpcnSigmaPos), static_cast<float>(cut.tpcnSigmaPos)}, "NSigma of positive particle axis binning in THnSparse."};
  ConfigurableAxis nsigmaaxisNeg{"nsigmaaxisNeg", {1, -static_cast<float>(cut.tpcnSigmaNeg), static_cast<float>(cut.tpcnSigmaNeg)}, "NSigma of negative particle axis binning in THnSparse."};

  // mixing
  using BinningTypeVzMu = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
  using BinningTypeVzCe = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  Configurable<int> numberofMixedEvents{"numberofMixedEvents", 5, "Number of events that should be mixed."};
  ConfigurableAxis axisVertexMixing{"axisVertexMixing", {5, -10, 10}, "Z vertex axis binning for mixing"};
  ConfigurableAxis axisMultiplicityMixing{"axisMultiplicityMixing", {5, 0, 5000}, "TPC multiplicity for bin"};
  ConfigurableAxis axisCentralityMixing{"axisCentralityMixing", {10, 0, 100}, "Multiplicity percentil binning for mixing"};

  // rotational
  Configurable<int> numberofRotations{"numberofRotations", 1, "Number of rotations for rotational background estimation."};
  Configurable<int> startingAngle{"startingAngle", 0, "Starting angle for rotational background estimation."};

  // other axes
  ConfigurableAxis axisNch{"axisNch", {1000, 0.0f, +1000.0f}, "Number of charged particles."};
  ConfigurableAxis axisResolutionPt{"axisResolutionPt", {1001, -1.0f, +1.0f}, "Resolution of Pt."};
  ConfigurableAxis axisResolutionPtPhi{"axisResolutionPtPhi", {1001, -0.01f, +0.01f}, "Resolution of Pt and Phi."};
  ConfigurableAxis axisResolutionMass{"axisResolutionMass", {1001, -0.01f, +0.01f}, "Resolution of Mass."};
  ConfigurableAxis axisResolutionVz{"axisResolutionVz", {1001, -3.0f, +3.0f}, "Resolution of Vz."};
  ConfigurableAxis massShiftAxis{"massShiftAxis", {1001, -0.02f, 0.02f}, "Mass correction axis."};

  // Axes specifications
  AxisSpec posZaxis = {400, -20., 20., "V_{z} (cm)"};
  AxisSpec dcaXYaxis = {1000, -1.0, 1.0, "DCA_{xy} (cm)"};
  AxisSpec dcaZaxis = {1000, -1.0, 1.0, "DCA_{z} (cm)"};
  AxisSpec etaQAaxis = {1000, -1.0, 1.0, "#eta"};
  AxisSpec tpcNClsFoundQAaxis = {110, 50., 160., "tpcNClsFound"};
  AxisSpec massShiftRelAxis = {101, -0.03f, 0.03f, ""};

  HistogramRegistry registry{"registry"};
  o2::analysis::rsn::Output* rsnOutput = nullptr;

  Service<o2::framework::O2DatabasePDG> pdg;

  int n = 0;
  float massPos = o2::track::PID::getMass(3);
  float massNeg = o2::track::PID::getMass(3);
  int pion = 2;
  int kaon = 3;
  int proton = 4;
  double* pointPair = nullptr;
  double* pointSys = nullptr;
  ROOT::Math::PxPyPzMVector d1, d2, mother, motherGen;
  bool produceTrue, produceLikesign, produceQA, produceStats, produceRotational, dataQA, MCTruthQA, globalTrack, inelGrater0, tpcPidOnly = false;
  float tpcnSigmaPos = 100.0f;
  float tpcnSigmaNeg = 100.0f;
  float combinedNSigma = 100.0f;
  float ptTOFThreshold = 0.5f;
  int tpcNClsFound = 70;
  int dauSize = 2;
  float vzCut = 10.0f;
  rsn::MixingType mixingType = rsn::MixingType::none;

  Filter triggerFilter = (o2::aod::evsel::sel8 == true);
  Filter vtxFilter = (nabs(o2::aod::collision::posZ) < vzCut);

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Ms>;
  using EventCandidate = EventCandidates::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr>;

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

    // Sparse axes
    AxisSpec invAxis = {invaxis, "Inv. mass (GeV/c^{2})", "im"};
    AxisSpec ptAxis = {ptaxis, "p_{T} (GeV/c)", "pt"};
    AxisSpec muAxis = {multiplicityaxis, "N", "mu"};
    AxisSpec mumAxis = {multiplicityaxis, "N", "mum"};
    AxisSpec ceAxis = {centralityaxis, "N", "ce"};
    AxisSpec cemAxis = {centralityaxis, "N", "cem"};
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

    produceQA = static_cast<bool>(produce.produceQA);
    produceStats = static_cast<bool>(produce.produceStats);
    produceTrue = static_cast<bool>(produce.produceTrue);
    produceLikesign = static_cast<bool>(produce.produceLikesign);
    mixingType = rsn::mixingTypeName(static_cast<std::string>(produce.eventMixing));
    produceRotational = static_cast<bool>(produce.produceRotational);
    tpcnSigmaPos = static_cast<float>(cut.tpcnSigmaPos);
    tpcnSigmaNeg = static_cast<float>(cut.tpcnSigmaNeg);
    tpcNClsFound = static_cast<int>(cut.tpcNClsFound);
    globalTrack = static_cast<bool>(cut.globalTrack);
    inelGrater0 = static_cast<bool>(cut.inelGrater0);
    combinedNSigma = static_cast<float>(cut.combinedNSigma);
    tpcPidOnly = static_cast<bool>(cut.tpcPidOnly);
    ptTOFThreshold = static_cast<float>(cut.ptTOFThreshold);
    vzCut = static_cast<float>(cut.vzCut);

    pointPair = new double[static_cast<int>(o2::analysis::rsn::PairAxisType::unknown)];
    pointSys = new double[static_cast<int>(o2::analysis::rsn::SystematicsAxisType::unknown)];
    rsnOutput = new o2::analysis::rsn::OutputSparse();
    rsnOutput->init(sparseAxes, allAxes, sysAxes, allAxesSys, produceTrue, mixingType, produceLikesign, produceRotational, &registry);

    // Print summary of configuration
    LOGF(info, "=== PhianalysisTHnSparse configuration summary ===");
    LOGF(info, "produceQA: %s", produceQA ? "true" : "false");
    LOGF(info, "produceStats: %s", produceStats ? "true" : "false");
    LOGF(info, "produceTrue: %s", static_cast<bool>(produce.produceTrue) ? "true" : "false");
    LOGF(info, "produceLikesign: %s", static_cast<bool>(produce.produceLikesign) ? "true" : "false");
    LOGF(info, "produceRotational: %s", static_cast<bool>(produce.produceRotational) ? "true" : "false");
    LOGF(info, "eventMixing: %s", static_cast<std::string>(produce.eventMixing).c_str());
    LOGF(info, "daughterPos: %d (PDG: %d)", static_cast<int>(daughterPos), static_cast<int>(daughterPosPDG));
    LOGF(info, "daughterNeg: %d (PDG: %d)", static_cast<int>(daughterNeg), static_cast<int>(daughterNegPDG));
    LOGF(info, "motherPDG: %d", static_cast<int>(motherPDG));
    LOGF(info, "tpcnSigmaPos: %.2f", tpcnSigmaPos);
    LOGF(info, "tpcnSigmaNeg: %.2f", tpcnSigmaNeg);
    LOGF(info, "tpcPidOnly: %s", tpcPidOnly ? "true" : "false");
    LOGF(info, "combinedNSigma: %.2f", combinedNSigma);
    LOGF(info, "ptTOFThreshold: %.2f", ptTOFThreshold);
    LOGF(info, "rapidity: %.2f", static_cast<float>(cut.rapidity));
    LOGF(info, "etatrack: %.2f", static_cast<float>(cut.etatrack));
    LOGF(info, "pt (min): %.2f", static_cast<float>(cut.pt));
    LOGF(info, "dcaXY: %.2f", static_cast<float>(cut.dcaXY));
    LOGF(info, "dcaZ: %.2f", static_cast<float>(cut.dcaZ));
    LOGF(info, "globalTrack: %s", globalTrack ? "true" : "false");
    LOGF(info, "inelGrater0: %s", inelGrater0 ? "true" : "false");
    LOGF(info, "tpcNClsFound: %d", tpcNClsFound);
    LOGF(info, "vzCut: %.2f", vzCut);
    LOGF(info, "mixingType: %d", static_cast<int>(mixingType));
    LOGF(info, "numberofMixedEvents: %d", static_cast<int>(numberofMixedEvents));
    LOGF(info, "numberofRotations: %d", static_cast<int>(numberofRotations));
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

    if (produceQA) {
      // Event QA
      registry.add("QAEvent/hSelection", "Event selection statistics", kTH1D, {{4, 0.0f, 4.0f}});
      auto hEvent = registry.get<TH1>(HIST("QAEvent/hSelection"));

      hEvent->GetXaxis()->SetBinLabel(1, "all events");
      hEvent->GetXaxis()->SetBinLabel(2, "Events passing trigger sel8");
      hEvent->GetXaxis()->SetBinLabel(3, "Events passing |V_{z}| cut");
      hEvent->GetXaxis()->SetBinLabel(4, "Events passing INEL>0 cut");
      hEvent->SetMinimum(0.1);

      registry.add("QAEvent/hVtxZ", "Vertex position along the z-axis", kTH1F, {posZaxis});
      registry.add("QAEvent/hCent", "Distribution of multiplicity percentile", kTH1F, {{101, 0., 101.}});
      registry.add("QAEvent/hMult", "Multiplicity (amplitude of non-zero channels in the FT0A + FT0C) ", kTH1F, {{300, 0., 30000.}});

      // Track QA
      registry.add("QATrack/hSelection", "Tracks statistics", kTH1D, {{9, 0.0f, 9.0f}});
      auto hTrack = registry.get<TH1>(HIST("QATrack/hSelection"));
      hTrack->GetXaxis()->SetBinLabel(1, "all tracks");
      hTrack->GetXaxis()->SetBinLabel(2, "passed pT cut");
      hTrack->GetXaxis()->SetBinLabel(3, "passed eta cut");
      hTrack->GetXaxis()->SetBinLabel(4, "passed DCA cut");
      hTrack->GetXaxis()->SetBinLabel(5, "passed PID cut");
      hTrack->GetXaxis()->SetBinLabel(6, "passed tpcNClsFound cut");
      hTrack->GetXaxis()->SetBinLabel(7, "passed isPrimaryTrack cut");
      hTrack->GetXaxis()->SetBinLabel(8, "passed isPVContributor cut");
      hTrack->GetXaxis()->SetBinLabel(9, "passed all cuts");
      hTrack->SetMinimum(0.1);

      registry.add("QATrack/hRapidity", "Rapidity distribution of K^{+} and K^{-}", kTH1F, {{200, -1, 1}});
      registry.add("QATrack/hEta", "Pseudorapidity distribution of K^{+} and K^{-}", kTH1F, {{200, -1, 1}});
      registry.add("QATrack/hTPCNClsFound", "Distribution of TPC NClsFound of K^{+} and K^{-}", kTH1F, {tpcNClsFoundQAaxis});
      registry.add("QATrack/hDCAxy", "Distribution of DCA_{xy} of K^{+} and K^{-}", kTH1F, {dcaXYaxis});
      registry.add("QATrack/hDCAz", "Distribution of DCA_{z} of K^{+} and K^{-}", kTH1F, {dcaZaxis});
      registry.add("QATrack/hPt", "Distribution of p_{T} of K^{+} and K^{-}", kTH1F, {ptaxis});

      // Phi candidate QA
      registry.add("QAPhi/hRapidity", "Rapidity distribution of #Phi candidates", kTH1F, {{200, -1, 1}});
      registry.add("QAPhi/hEta", "Pseudorapidity distribution of #Phi candidates", kTH1F, {{200, -1, 1}});
      registry.add("QAPhi/hdPhi", "Azimuthal distribution (#Delta#phi) of #Phi candidates", kTH1F, {{100, -o2::constants::math::TwoPI, o2::constants::math::TwoPI}});
      auto hdPhi = registry.get<TH1>(HIST("QAPhi/hdPhi"));
      hdPhi->GetXaxis()->SetTitle("#Delta#phi (rad)");

      registry.add("QAPhi/h2dPhiPt", "Azimuthal distribution (#Delta#phi) of #Phi candidates vs p_{T}", kTH2F, {ptaxis, {100, -o2::constants::math::TwoPI, o2::constants::math::TwoPI}});
      auto h2dPhiPt = registry.get<TH2>(HIST("QAPhi/h2dPhiPt"));
      h2dPhiPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h2dPhiPt->GetYaxis()->SetTitle("#Delta#phi (rad)");

      registry.add("QAPhi/hTheta", "Polar distribution of #Phi candidates", kTH1F, {{100, 0.0f, o2::constants::math::PI}});
      auto hTheta = registry.get<TH1>(HIST("QAPhi/hTheta"));
      hTheta->GetXaxis()->SetTitle("#theta (rad)");

      registry.add("QAPhi/h2dThetaPt", "Polar distribution (#Delta#theta) of #Phi candidates vs p_{T}", kTH2F, {ptaxis, {100, -o2::constants::math::PI, o2::constants::math::PI}});

      auto h2dThetaPt = registry.get<TH2>(HIST("QAPhi/h2dThetaPt"));
      h2dThetaPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h2dThetaPt->GetYaxis()->SetTitle("#Delta#theta (rad)");

      // Rotational background QA
      if (produceRotational) {
        registry.add("QARotational/hRapidity", "Rapidity distribution of #Phi candidates from rotational background", kTH1F, {{200, -1, 1}});
        registry.add("QARotational/hEta", "Pseudorapidity distribution of #Phi candidates from rotational background", kTH1F, {{200, -1, 1}});
        registry.add("QARotational/hdPhi", "Rotational background: Azimuthal distribution (#Delta#phi)", kTH1F, {{100, -o2::constants::math::TwoPI, o2::constants::math::TwoPI}});
        auto hRPhi = registry.get<TH1>(HIST("QARotational/hdPhi"));
        hRPhi->GetXaxis()->SetTitle("#Delta#phi");

        registry.add("QARotational/h2dPhiPt", "Rotational background: Azimuthal distribution (#Delta#phi) vs p_{T}", kTH2F, {ptaxis, {100, -o2::constants::math::TwoPI, o2::constants::math::TwoPI}});
        auto hR2dPhiPt = registry.get<TH2>(HIST("QARotational/h2dPhiPt"));
        hR2dPhiPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hR2dPhiPt->GetYaxis()->SetTitle("#Delta#phi");

        registry.add("QARotational/hTheta", "Rotational background: Polar distribution (#theta)", kTH1F, {{100, 0.0f, o2::constants::math::PI}});
        auto hRdTheta = registry.get<TH1>(HIST("QARotational/hTheta"));
        hRdTheta->GetXaxis()->SetTitle("#theta (rad)");

        registry.add("QARotational/h2dThetaPt", "Rotational background: Polar distribution (#Delta#theta) vs p_{T}", kTH2F, {ptaxis, {100, -o2::constants::math::PI, o2::constants::math::PI}});
        auto hR2dThetaPt = registry.get<TH2>(HIST("QARotational/h2dThetaPt"));
        hR2dThetaPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hR2dThetaPt->GetYaxis()->SetTitle("#Delta#theta");
      }

      // Mixing QA
      if (mixingType != rsn::MixingType::none) {
        registry.add("QAMixing/hSelection", "Event mixing selection statistics", kTH1D, {{1, 0.0f, 1.0f}});
        auto hEM = registry.get<TH1>(HIST("QAMixing/hSelection"));
        hEM->GetXaxis()->SetBinLabel(1, "Full event mixing statistics");
        hEM->SetMinimum(0.1);

        registry.add("QAMixing/h2mu1_mu2", "Event Mixing Multiplicity", kTH2F, {axisMultiplicityMixing, axisMultiplicityMixing});
        auto h2EMmu = registry.get<TH2>(HIST("QAMixing/h2mu1_mu2"));
        h2EMmu->GetXaxis()->SetTitle("1.Event multiplicity");
        h2EMmu->GetYaxis()->SetTitle("2.Event multiplicity");

        registry.add("QAMixing/h2ce1_ce2", "Event Mixing Centrality", kTH2F, {axisCentralityMixing, axisCentralityMixing});
        auto h2EMce = registry.get<TH2>(HIST("QAMixing/h2ce1_ce2"));
        h2EMce->GetXaxis()->SetTitle("1.Event centrality");
        h2EMce->GetYaxis()->SetTitle("2.Event centrality");

        registry.add("QAMixing/h2vz1_vz2", "Event Mixing Vertex z", kTH2F, {axisVertexMixing, axisVertexMixing});
        auto hEMTvz = registry.get<TH2>(HIST("QAMixing/h2vz1_vz2"));
        hEMTvz->GetXaxis()->SetTitle("1.Event V_{z}");
        hEMTvz->GetYaxis()->SetTitle("2.Event V_{z}");
      }

      // PID QA
      // TPC
      registry.add("QAPID/hTPCnSigma", "Distribution of TPC nSigma of K^{+} and K^{-}", kTH1F, {{200, -10, 10}});
      auto hTPCnSigma = registry.get<TH1>(HIST("QAPID/hTPCnSigma"));
      hTPCnSigma->GetXaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");

      registry.add("QAPID/h2TPCnSigma", "", kTH2F, {{200, -10, 10}, {200, -10, 10}});
      auto h2TPCnSigma = registry.get<TH2>(HIST("QAPID/h2TPCnSigma"));
      h2TPCnSigma->GetXaxis()->SetTitle("n#sigma_{TPC} K^{+}");
      h2TPCnSigma->GetYaxis()->SetTitle("n#sigma_{TPC} K^{-}");

      registry.add("QAPID/h2TPCnSigmaPt", "", kTH2F, {ptaxis, {200, -10, 10}});
      auto h2TPCnSigmaPt = registry.get<TH2>(HIST("QAPID/h2TPCnSigmaPt"));
      h2TPCnSigmaPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h2TPCnSigmaPt->GetYaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");

      // TOF
      registry.add("QAPID/hTOFnSigma", "Distribution of TOF nSigma of K^{+} and K^{-}", kTH1F, {{200, -10, 10}});
      auto hTOFnSigma = registry.get<TH1>(HIST("QAPID/hTOFnSigma"));
      hTOFnSigma->GetXaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");

      registry.add("QAPID/h2TOFnSigma", "", kTH2F, {{200, -10, 10}, {200, -10, 10}});
      auto h2TOFnSigma = registry.get<TH2>(HIST("QAPID/h2TOFnSigma"));
      h2TOFnSigma->GetXaxis()->SetTitle("n#sigma_{TOF} K^{+}");
      h2TOFnSigma->GetYaxis()->SetTitle("n#sigma_{TOF} K^{-}");

      registry.add("QAPID/h2TOFnSigmaPt", "", kTH2F, {ptaxis, {200, -10, 10}});
      auto h2TOFnSigmaPt = registry.get<TH2>(HIST("QAPID/h2TOFnSigmaPt"));
      h2TOFnSigmaPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h2TOFnSigmaPt->GetYaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");

      // MC
      if (static_cast<bool>(produce.produceTrue)) {
        // Rec
        registry.add("QAMC/Rec/hSelection", "MC Rec Event statistics", kTH1F, {{5, 0.0f, 5.0f}});
        auto hMCEventTruth = registry.get<TH1>(HIST("QAMC/Rec/hSelection"));
        hMCEventTruth->GetXaxis()->SetBinLabel(1, "Full MC Rec event statistics");
        hMCEventTruth->GetXaxis()->SetBinLabel(2, "MC Rec events passing sel8 cut");
        hMCEventTruth->GetXaxis()->SetBinLabel(3, "MC Rec events passing V_{z} cut");
        hMCEventTruth->GetXaxis()->SetBinLabel(4, "MC Rec events with V_{z} cut and INEL>0");
        hMCEventTruth->GetXaxis()->SetBinLabel(5, "Reconstructed #Phi candidates matched to true #Phi");
        hMCEventTruth->SetMinimum(0.1);

        registry.add("QAMC/Rec/hInvMassTrueFalse", "", kTH1F, {invAxis}); // not written events in True distribution due to repetition of mothers

        // Gen
        registry.add("QAMC/Gen/hSelection", "MC Gen Event statistics", kTH1F, {{4, 0.0f, 4.0f}});
        auto hMCEventGen = registry.get<TH1>(HIST("QAMC/Gen/hSelection"));
        hMCEventGen->GetXaxis()->SetBinLabel(1, "Full MC Gen event statistics");
        hMCEventGen->GetXaxis()->SetBinLabel(2, "MC Gen events within V_{z} cut");
        hMCEventGen->GetXaxis()->SetBinLabel(3, "MC Gen events with V_{z} cut and INEL>0");
        hMCEventGen->GetXaxis()->SetBinLabel(4, "Generated #Phi candidates");
        hMCEventGen->SetMinimum(0.1);

        // Resolution
        registry.add("Factors/h2ResolutionVz", "Resolution of collision V_{z}", kTH2F, {vzaxis, axisResolutionVz});
        auto hResVz = registry.get<TH2>(HIST("Factors/h2ResolutionVz"));
        hResVz->GetXaxis()->SetTitle("V_{z}^{rec} (cm)");
        hResVz->GetYaxis()->SetTitle("#DeltaV_{z} = V_{z}^{rec} - V_{z}^{gen} (cm)");
        registry.add("Factors/h2ResolutionPt", "Resolution of charged particles p_{T}", kTH2F, {ptaxis, axisResolutionPt});
        auto hResPt = registry.get<TH2>(HIST("Factors/h2ResolutionPt"));

        hResPt->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
        hResPt->GetYaxis()->SetTitle("#Deltap_{T} = p_{T}^{rec} - p_{T}^{gen} (GeV/c)");
        registry.add("Factors/h2ResolutionPtPhi", "p_{T} resolution vs p_{T}^{rec}", kTH2F, {ptaxis, axisResolutionPtPhi});
        auto hResPtPhi = registry.get<TH2>(HIST("Factors/h2ResolutionPtPhi"));
        hResPtPhi->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
        hResPtPhi->GetYaxis()->SetTitle("#Deltap_{T} = p_{T}^{rec} - p_{T}^{gen} (GeV/c)");

        registry.add("Factors/h2MassResolution", "Mass resolution vs p_{T}^{rec}", kTH2F, {ptaxis, axisResolutionMass});
        auto hResMass = registry.get<TH2>(HIST("Factors/h2MassResolution"));
        hResMass->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
        hResMass->GetYaxis()->SetTitle("#Deltam = m^{gen}_{KK} - m^{rec}_{KK} (GeV/c^{2})");

        registry.add("Factors/h2MassShift", "Mass shift vs p_{T}^{rec}", kTH2F, {ptaxis, massShiftAxis});
        auto hResMassGen = registry.get<TH2>(HIST("Factors/h2MassShift"));
        hResMassGen->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
        hResMassGen->GetYaxis()->SetTitle("#Deltam = m^{gen}_{#phi} - m^{gen}_{KK} (GeV/c^{2})");

        registry.add("Factors/h2MassShiftRel", "Relative mass shift vs p_{T}^{rec}", kTH2F, {ptaxis, massShiftRelAxis});
        auto hMassCorr = registry.get<TH2>(HIST("Factors/h2MassShiftRel"));
        hMassCorr->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
        hMassCorr->GetYaxis()->SetTitle("m^{gen}_{#phi} - m^{gen}_{KK}/m^{gen}_{#phi}");
      }
    }
    // Factors
    registry.add("Factors/hCentralityVsMultMC", "Event centrality vs MC multiplicity", kTH2F, {{101, 0.0f, 101.0f}, axisNch});
    registry.add("Factors/hCentralityVsMult", "Event centrality vs multiplicity", kTH2F, {{101, 0.0f, 101.0f}, axisNch});
    registry.add("Factors/hEventCentrality", "Event centrality", kTH1F, {{101, 0, 101}});
    registry.add("Factors/hNrecInGen", "Number of collisions in MC", kTH1F, {{3, -0.5, 2.5}});
    registry.add("Factors/hGenEvents", "Generated events", HistType::kTH2F, {{axisNch}, {4, 0, 4}});
    auto hGenEvents = registry.get<TH2>(HIST("Factors/hGenEvents"));
    hGenEvents->GetYaxis()->SetBinLabel(1, "All generated events");
    hGenEvents->GetYaxis()->SetBinLabel(2, "Generated events passing V_{z} cut");
    hGenEvents->GetYaxis()->SetBinLabel(3, "Generated events passing INEL>0");
    hGenEvents->GetYaxis()->SetBinLabel(4, "Generated events with at least one reconstructed event");
    registry.add("Factors/h2dGenPhi", "Centrality vs p_{T}", kTH2D, {{101, 0.0f, 101.0f}, ptaxis});
    registry.add("Factors/h3dGenPhiVsMultMCVsCentrality", "MC multiplicity vs centrality vs p_{T}", kTH3D, {axisNch, {101, 0.0f, 101.0f}, ptaxis});
  }

  template <typename T>
  float tpcNsigma(const T& track)
  {
    float tpcNsigma = 0.0f;
    int particleType = (track.sign() > 0) ? static_cast<int>(daughterPos) : static_cast<int>(daughterNeg);

    if (particleType == pion)
      tpcNsigma = track.tpcNSigmaPi();
    else if (particleType == kaon)
      tpcNsigma = track.tpcNSigmaKa();
    else if (particleType == proton)
      tpcNsigma = track.tpcNSigmaPr();
    return tpcNsigma;
  }
  template <typename T>
  float tofNsigma(const T& track)
  {
    float tofNsigma = 0.0f;
    int particleType = (track.sign() > 0) ? static_cast<int>(daughterPos) : static_cast<int>(daughterNeg);

    if (particleType == pion)
      tofNsigma = track.tofNSigmaPi();
    else if (particleType == kaon)
      tofNsigma = track.tofNSigmaKa();
    else if (particleType == proton)
      tofNsigma = track.tofNSigmaPr();
    return tofNsigma;
  }
  template <typename T>
  bool selectedTrack(const T& track, bool isPositive)
  {
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 0.5); // all tracks

    // Apply pT cut
    if (track.pt() < static_cast<float>(cut.pt))
      return false;
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 1.5);

    // Apply eta cut
    if (std::abs(track.eta()) >= static_cast<float>(cut.etatrack))
      return false;
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 2.5);

    // Apply DCA cuts
    if (std::abs(track.dcaXY()) >= static_cast<float>(cut.dcaXY) ||
        std::abs(track.dcaZ()) >= static_cast<float>(cut.dcaZ))
      return false;
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 3.5);

    // PID selection: TPC-only for pt < threshold value, TPC+TOF for pt >= threshold value and have TOF, else TPC-only
    float nSigmaCut = isPositive ? tpcnSigmaPos : tpcnSigmaNeg;
    if (track.pt() < ptTOFThreshold || !track.hasTOF() || tpcPidOnly) {
      if (std::abs(tpcNsigma(track)) >= nSigmaCut)
        return false;
    } else {
      if (std::sqrt(tpcNsigma(track) * tpcNsigma(track) + tofNsigma(track) * tofNsigma(track)) >= combinedNSigma)
        return false;
    }
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 4.5);

    // Apply tpcNClsFound cut
    if (track.tpcNClsFound() < tpcNClsFound)
      return false;
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 5.5);

    if (globalTrack) {
      // Apply Global track cuts
      if (!track.isGlobalTrack())
        return false;
    } else {
      // Apply Primary track cuts
      if (!track.isPrimaryTrack())
        return false;
    }
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 6.5);

    // Apply PV Contributor cuts
    if (!track.isPVContributor())
      return false;
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 7.5);

    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 8.5);

    return true;
  }
  template <typename T>
  ROOT::Math::PxPyPzMVector calculateMother(const T& track1, const T& track2)
  {
    d1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPos);
    d2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massNeg);
    return d1 + d2;
  }
  bool seletectedMother(const ROOT::Math::PxPyPzMVector& mother)
  {
    if (std::abs(mother.Rapidity()) > static_cast<float>(cut.rapidity))
      return false;
    return true;
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

  void processData(EventCandidate const& collision, TrackCandidates const& /*tracks*/)
  {
    auto posDaughters = positive->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDaughters = negative->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    int nch = 0;

    if (produceQA)
      registry.fill(HIST("QAEvent/hSelection"), 0.5);

    if (!collision.sel8())
      return;

    if (produceQA)
      registry.fill(HIST("QAEvent/hSelection"), 1.5);

    if (std::abs(collision.posZ()) > vzCut)
      return;

    if (produceQA)
      registry.fill(HIST("QAEvent/hSelection"), 2.5);

    if (inelGrater0 && !collision.isInelGt0())
      return;

    registry.fill(HIST("Factors/hEventCentrality"), collision.centFT0M());

    if (produceQA) {
      registry.fill(HIST("QAEvent/hSelection"), 3.5);
      registry.fill(HIST("QAEvent/hVtxZ"), collision.posZ());
      registry.fill(HIST("QAEvent/hMult"), getMultiplicity(collision));
      registry.fill(HIST("QAEvent/hCent"), getCentrality(collision));

      if (produceStats) {
        dataQA = true;
        for (const auto& track : posDaughters) {
          if (track.isPrimaryTrack() && std::abs(track.eta()) >= static_cast<float>(cut.etatrack))
            nch++;
          selectedTrack(track, true);
        }
        for (const auto& track : negDaughters) {
          if (track.isPrimaryTrack() && std::abs(track.eta()) >= static_cast<float>(cut.etatrack))
            nch++;
          selectedTrack(track, false);
        }
        dataQA = false;
      }
      registry.fill(HIST("Factors/hCentralityVsMult"), getCentrality(collision), nch);
    }

    if (static_cast<int>(verbose.verboselevel) > 0 && static_cast<int>(verbose.refresh) > 0 && collision.globalIndex() % static_cast<int>(verbose.refresh) == static_cast<int>(verbose.refreshIndex))
      LOGF(info, "%d pos=%lld neg=%lld, Z vertex position: %f [cm]", collision.globalIndex(), posDaughters.size(), negDaughters.size(), collision.posZ());

    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughters, negDaughters))) {

      if (!selectedTrack(track1, true)) // track1 is positive
        continue;
      if (!selectedTrack(track2, false)) // track2 is negative
        continue;

      mother = calculateMother(track1, track2);
      if (!seletectedMother(mother))
        continue;

      if (produceQA) {
        registry.fill(HIST("QAPID/h2TPCnSigma"), tpcNsigma(track1), tpcNsigma(track2));
        registry.fill(HIST("QAPID/h2TPCnSigmaPt"), track1.pt(), tpcNsigma(track1));
        registry.fill(HIST("QAPID/h2TPCnSigmaPt"), track2.pt(), tpcNsigma(track2));

        registry.fill(HIST("QAPID/h2TOFnSigma"), tofNsigma(track1), tofNsigma(track2));
        registry.fill(HIST("QAPID/h2TOFnSigmaPt"), track1.pt(), tofNsigma(track1));
        registry.fill(HIST("QAPID/h2TOFnSigmaPt"), track2.pt(), tofNsigma(track2));

        registry.fill(HIST("QAPID/hTPCnSigma"), tpcNsigma(track1));
        registry.fill(HIST("QAPID/hTPCnSigma"), tpcNsigma(track2));
        if (track1.hasTOF())
          registry.fill(HIST("QAPID/hTOFnSigma"), tofNsigma(track1));
        if (track2.hasTOF())
          registry.fill(HIST("QAPID/hTOFnSigma"), tofNsigma(track2));

        registry.fill(HIST("QATrack/hEta"), track1.eta());
        registry.fill(HIST("QATrack/hEta"), track2.eta());
        registry.fill(HIST("QATrack/hPt"), track1.pt());
        registry.fill(HIST("QATrack/hPt"), track2.pt());
        registry.fill(HIST("QATrack/hDCAxy"), track1.dcaXY());
        registry.fill(HIST("QATrack/hDCAxy"), track2.dcaXY());
        registry.fill(HIST("QATrack/hDCAz"), track1.dcaZ());
        registry.fill(HIST("QATrack/hDCAz"), track2.dcaZ());
        registry.fill(HIST("QATrack/hTPCNClsFound"), track1.tpcNClsFound());
        registry.fill(HIST("QATrack/hTPCNClsFound"), track2.tpcNClsFound());
        registry.fill(HIST("QATrack/hRapidity"), track1.rapidity(massPos));
        registry.fill(HIST("QATrack/hRapidity"), track2.rapidity(massNeg));

        registry.fill(HIST("QAPhi/hRapidity"), mother.Rapidity());
        registry.fill(HIST("QAPhi/hEta"), mother.Eta());
        registry.fill(HIST("QAPhi/hdPhi"), track1.phi() - track2.phi());
        registry.fill(HIST("QAPhi/h2dPhiPt"), mother.Pt(), track1.phi() - track2.phi());
        registry.fill(HIST("QAPhi/hTheta"), mother.Theta());
        registry.fill(HIST("QAPhi/h2dThetaPt"), mother.Pt(), d1.Theta() - d2.Theta());
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
      rsnOutput->fillUnlikepm(pointPair);

      if (produceRotational) {
        for (int i = 1; i <= static_cast<int>(numberofRotations); i++) {
          float starting = static_cast<float>(startingAngle) * TMath::DegToRad();
          float angle = starting + i * ((o2::constants::math::TwoPI - 2 * starting) / (static_cast<int>(numberofRotations) + 1));
          float px2new = track2.px() * std::cos(angle) - track2.py() * std::sin(angle);
          float py2new = track2.px() * std::sin(angle) + track2.py() * std::cos(angle);
          d2 = ROOT::Math::PxPyPzMVector(px2new, py2new, track2.pz(), massNeg);
          mother = d1 + d2;

          if (produceQA) {
            registry.fill(HIST("QARotational/hRapidity"), mother.Rapidity());
            registry.fill(HIST("QARotational/hEta"), mother.Eta());
            registry.fill(HIST("QARotational/hdPhi"), d1.Phi() - d2.Phi());
            registry.fill(HIST("QARotational/h2dPhiPt"), mother.Pt(), d1.Phi() - d2.Phi());
            registry.fill(HIST("QARotational/hTheta"), mother.Theta());
            registry.fill(HIST("QARotational/h2dThetaPt"), mother.Pt(), d1.Theta() - d2.Theta());
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

          rsnOutput->fillRotationpm(pointPair);
        }
      }
    }

    if (produceLikesign) {

      for (const auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posDaughters, posDaughters))) {
        if (!selectedTrack(track1, true)) // both positive
          continue;
        if (!selectedTrack(track2, true)) // both positive
          continue;

        mother = calculateMother(track1, track2);
        if (!seletectedMother(mother))
          continue;

        if (static_cast<int>(verbose.verboselevel) > 1)
          LOGF(info, "Like-sign positive: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.M());

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
        if (!selectedTrack(track1, false)) // both negative
          continue;
        if (!selectedTrack(track2, false)) // both negative
          continue;

        mother = calculateMother(track1, track2);
        if (!seletectedMother(mother))
          continue;

        if (static_cast<int>(verbose.verboselevel) > 1)
          LOGF(info, "Like-sign negative: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.M());

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
    if (!static_cast<bool>(produce.produceTrue))
      return;

    registry.fill(HIST("QAMC/Rec/hSelection"), 0.5);

    if (!collision.sel8())
      return;

    if (produceQA)
      registry.fill(HIST("QAMC/Rec/hSelection"), 1.5);

    auto posDaughtersMC = positiveMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDaughtersMC = negativeMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (!collision.has_mcCollision()) {
      if (static_cast<int>(verbose.verboselevel) > 0)
        LOGF(warning, "No MC collision for this collision, skip...");
      return;
    }
    auto mcCollision = collision.mcCollision();
    registry.fill(HIST("Factors/h2ResolutionVz"), collision.posZ(), (collision.posZ() - mcCollision.posZ()));

    if (std::abs(mcCollision.posZ()) > vzCut)
      return;

    if (produceQA)
      registry.fill(HIST("QAMC/Rec/hSelection"), 2.5);

    if (inelGrater0 && !collision.isInelGt0())
      return;

    if (produceQA)
      registry.fill(HIST("QAMC/Rec/hSelection"), 3.5);

    for (const auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto mctrack = track.mcParticle();
        registry.fill(HIST("Factors/h2ResolutionPt"), track.pt(), (track.pt() - mctrack.pt()));
      }
    }

    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersMC, negDaughtersMC))) {

      if (!track1.has_mcParticle()) {
        if (static_cast<int>(verbose.verboselevel) > 0)
          LOGF(warning, "No MC particle for track, skip...");
        continue;
      }

      if (!track2.has_mcParticle()) {
        if (static_cast<int>(verbose.verboselevel) > 0)
          LOGF(warning, "No MC particle for track, skip...");
        continue;
      }

      if (!selectedTrack(track1, true)) // track1 is positive
        continue;
      if (!selectedTrack(track2, false)) // track2 is negative
        continue;

      const auto mctrack1 = track1.mcParticle();
      const auto mctrack2 = track2.mcParticle();
      int track1PDG = std::abs(mctrack1.pdgCode());
      int track2PDG = std::abs(mctrack2.pdgCode());

      if (!(track1PDG == daughterPosPDG && track2PDG == daughterNegPDG)) {
        continue;
      }
      n = 0;
      for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
        for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {

          if (mothertrack1.pdgCode() != mothertrack2.pdgCode())
            continue;

          if (mothertrack1.globalIndex() != mothertrack2.globalIndex())
            continue;

          if (std::abs(mothertrack1.y()) > static_cast<float>(cut.rapidity))
            continue;

          if (std::abs(mothertrack1.pdgCode()) != motherPDG)
            continue;

          mother = calculateMother(track1, track2);
          motherGen = calculateMother(mctrack1, mctrack2);
          if (!seletectedMother(mother))
            continue;

          if (n > 0) {
            if (produceQA)
              registry.fill(HIST("QAMC/Rec/hInvMassTrueFalse"), mother.M());
            continue;
          }

          if (static_cast<int>(verbose.verboselevel) > 1) {
            LOGF(info, "Collision: %ld True: %d, d1=%d (%ld), d2=%d (%ld), mother=%d (%ld)", collision.globalIndex(), n, mctrack1.pdgCode(), mctrack1.globalIndex(), mctrack2.pdgCode(), mctrack2.globalIndex(), mothertrack1.pdgCode(), mothertrack1.globalIndex());
            LOGF(info, "Track %d px: %f, py=%f, pz=%f, px: %f, py=%f, pz=%f", n, track1.px(), track1.py(), track1.pz(), track2.px(), track2.py(), track2.pz());
            LOGF(info, "mcTrack %d px: %f, py=%f, pz=%f, px: %f, py=%f, pz=%f", n, mctrack1.px(), mctrack1.py(), mctrack1.pz(), mctrack2.px(), mctrack2.py(), mctrack2.pz());
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

          if (produceQA)
            registry.fill(HIST("QAMC/Rec/hSelection"), 4.5);

          auto phiP = mothertrack1.p();
          auto phiE = mothertrack1.e();
          auto massGen = std::sqrt(phiE * phiE - phiP * phiP);

          registry.fill(HIST("Factors/h2ResolutionPtPhi"), mother.Pt(), (mother.Pt() - mothertrack1.pt()));
          registry.fill(HIST("Factors/h2MassResolution"), mother.Pt(), (motherGen.M() - mother.M()));
          registry.fill(HIST("Factors/h2MassShift"), mother.Pt(), (massGen - motherGen.M()));
          registry.fill(HIST("Factors/h2MassShiftRel"), mother.Pt(), (massGen - motherGen.M()) / massGen);

          if (static_cast<int>(verbose.verboselevel) > 1)
            LOGF(info, "mother.M()=%f, motherGen.M()=%f, massGen =%f", mother.M(), motherGen.M(), massGen);

          rsnOutput->fillUnliketrue(pointPair);

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

          rsnOutput->fillUnlikegenOld(pointPair);

          n++;
        }
      }
    }
  }
  PROCESS_SWITCH(PhianalysisTHnSparse, processTrue, "Process Event for MC reconstruction.", false);

  void processGen(McCollisionMults::iterator const& mcCollision, soa::SmallGroups<EventCandidatesMCGen> const& collisions, LabeledTracks const& /*particles*/, aod::McParticles const& mcParticles)
  {
    if (!static_cast<bool>(produce.produceTrue))
      return;

    if (produceQA)
      registry.fill(HIST("QAMC/Gen/hSelection"), 0.5);

    if (std::abs(mcCollision.posZ()) > vzCut)
      return;

    if (produceQA)
      registry.fill(HIST("QAMC/Gen/hSelection"), 1.5);

    if (inelGrater0 && !mcCollision.isInelGt0())
      return;

    if (produceQA)
      registry.fill(HIST("QAMC/Gen/hSelection"), 2.5);

    if (collisions.size() == 0)
      return;

    for (const auto& collision : collisions) {
      auto centralityGen = getCentrality(collision);
      auto multiplicityGen = getMultiplicity(collision);

      for (const auto& particle : mcParticles) {

        if (std::abs(particle.y()) > static_cast<float>(cut.rapidity))
          continue;

        if (particle.pdgCode() == motherPDG) {

          auto daughters = particle.daughters_as<aod::McParticles>();
          if (daughters.size() != dauSize)
            continue;

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
          if (!daup || !daun)
            continue;

          mother = d1 + d2;

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    multiplicityGen,
                                    centralityGen,
                                    0,
                                    0,
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    mcCollision.posZ(),
                                    0,
                                    0,
                                    0);

          rsnOutput->fillUnlikegen(pointPair);
          if (produceQA)
            registry.fill(HIST("QAMC/Gen/hSelection"), 3.5);
        }
      }
    }
  }
  PROCESS_SWITCH(PhianalysisTHnSparse, processGen, "Process MC Generated.", false);

  void processMixed(soa::Filtered<EventCandidates> const& collisions, TrackCandidates const& tracks)
  {
    if (mixingType == rsn::MixingType::none)
      return;

    auto tracksTuple = std::make_tuple(tracks);

    BinningTypeVzCe binningVzCe{{axisVertexMixing, axisCentralityMixing}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVzCe> pairVzCe{binningVzCe, static_cast<int>(numberofMixedEvents), -1, collisions, tracksTuple, &cache};

    BinningTypeVzMu binningVzMu{{axisVertexMixing, axisMultiplicityMixing}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVzMu> pairVzMu{binningVzMu, static_cast<int>(numberofMixedEvents), -1, collisions, tracksTuple, &cache};

    if (mixingType == rsn::MixingType::ce) {
      for (const auto& [c1, tracks1, c2, tracks2] : pairVzCe) {
        if (produceQA)
          registry.fill(HIST("QAMixing/hSelection"), 0.5);

        auto posDaughtersc1 = positive->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto posDaughtersc2 = positive->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);
        auto negDaughtersc1 = negative->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto negDaughtersc2 = negative->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);

        if (produceQA) {
          registry.fill(HIST("QAMixing/h2mu1_mu2"), getMultiplicity(c1), getMultiplicity(c2));
          registry.fill(HIST("QAMixing/h2ce1_ce2"), getCentrality(c1), getCentrality(c2));
          registry.fill(HIST("QAMixing/h2vz1_vz2"), c1.posZ(), c2.posZ());
        }

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersc1, negDaughtersc2))) {

          if (!selectedTrack(track1, true)) // track1 is positive
            continue;
          if (!selectedTrack(track2, false)) // track2 is negative
            continue;

          mother = calculateMother(track1, track2);
          if (!seletectedMother(mother))
            continue;

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

          if (!selectedTrack(track1, true)) // track1 is positive
            continue;
          if (!selectedTrack(track2, false)) // track2 is negative
            continue;

          mother = calculateMother(track1, track2);
          if (!seletectedMother(mother))
            continue;

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
        if (produceQA)
          registry.fill(HIST("QAMixing/hSelection"), 0.5);

        auto posDaughtersc1 = positive->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto posDaughtersc2 = positive->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);
        auto negDaughtersc1 = negative->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto negDaughtersc2 = negative->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);

        if (produceQA) {
          registry.fill(HIST("QAMixing/h2mu1_mu2"), getMultiplicity(c1), getMultiplicity(c2));
          registry.fill(HIST("QAMixing/h2ce1_ce2"), getCentrality(c1), getCentrality(c2));
          registry.fill(HIST("QAMixing/h2vz1_vz2"), c1.posZ(), c2.posZ());
        }

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersc1, negDaughtersc2))) {

          if (!selectedTrack(track1, true)) // track1 is positive
            continue;

          if (!selectedTrack(track2, false)) // track2 is negative
            continue;

          mother = calculateMother(track1, track2);
          if (!seletectedMother(mother))
            continue;

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

          if (!selectedTrack(track1, true))

            continue;
          if (!selectedTrack(track2, false))
            continue;

          mother = calculateMother(track1, track2);
          if (!seletectedMother(mother))
            continue;

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

  void processFactors(McCollisionMults::iterator const& mcCollision, soa::SmallGroups<EventCandidatesMCGen> const& collisions, LabeledTracks const& /*particles*/, aod::McParticles const& mcParticles)
  {
    registry.fill(HIST("Factors/hGenEvents"), mcCollision.multMCNParticlesEta08(), 0.5);

    if (std::abs(mcCollision.posZ()) > vzCut)
      return;

    registry.fill(HIST("Factors/hGenEvents"), mcCollision.multMCNParticlesEta08(), 1.5);

    if (inelGrater0 && !mcCollision.isInelGt0())
      return;

    registry.fill(HIST("Factors/hGenEvents"), mcCollision.multMCNParticlesEta08(), 2.5);

    float centrality = 100.5f;
    for (auto const& collision : collisions) {
      centrality = collision.centFT0M();
    }

    registry.fill(HIST("Factors/hCentralityVsMultMC"), centrality, mcCollision.multMCNParticlesEta08());
    registry.fill(HIST("Factors/hNrecInGen"), collisions.size());

    for (const auto& particle : mcParticles) {

      if (std::abs(particle.y()) > static_cast<float>(cut.rapidity))
        continue;

      if (particle.pdgCode() == motherPDG) {

        auto daughters = particle.daughters_as<aod::McParticles>();
        if (daughters.size() != dauSize)
          continue;

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
        if (!daup || !daun)
          continue;

        mother = d1 + d2;

        registry.fill(HIST("Factors/h2dGenPhi"), centrality, mother.Pt());
        registry.fill(HIST("Factors/h3dGenPhiVsMultMCVsCentrality"), mcCollision.multMCNParticlesEta08(), centrality, mother.Pt());
      }
    }

    if (collisions.size() == 0)
      return;

    registry.fill(HIST("Factors/hGenEvents"), mcCollision.multMCNParticlesEta08(), 3.5);
  }
  PROCESS_SWITCH(PhianalysisTHnSparse, processFactors, "Process to obtain normalization factors from MC.", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PhianalysisTHnSparse>(cfgc)};
}
