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
#include "Common/DataModel/PIDResponse.h"
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
    Configurable<bool> produceTrue{"produceTrue", false, "Produce True and Gen histograms."};
    Configurable<bool> produceLikesign{"produceLikesign", false, "Produce Like sign histograms."};
    Configurable<std::string> eventMixing{"eventMixing", "none", "Produce Event Mixing histograms of type."};
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
    Configurable<float> etapair{"etapair", 0.5f, "Eta cut for pair."};
    Configurable<float> pt{"pt", 0.15f, "Cut: Minimal value of tracks pt."};
    Configurable<float> dcaXY{"dcaXY", 1.0f, "Cut: Maximal value of tracks DCA XY."};
    Configurable<float> dcaZ{"dcaZ", 1.0f, "Cut: Maximal value of tracks DCA Z."};
    Configurable<bool> globalTrack{"globalTrack", false, "Use global track selection."};
    Configurable<int> tpcNClsFound{"tpcNClsFound", 70, "Cut: Minimal value of found TPC clasters"};
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

  // Axes specifications
  AxisSpec posZaxis = {400, -20., 20., "V_{z} (cm)"};
  AxisSpec dcaXYaxis = {800, -2.0, 2.0, "DCA_{xy} (cm)"};
  AxisSpec dcaZaxis = {800, -2.0, 2.0, "DCA_{z} (cm)"};
  AxisSpec etaQAaxis = {1000, -1.0, 1.0, "#eta"};
  AxisSpec tpcNClsFoundQAaxis = {110, 50., 160., "tpcNClsFound"};

  HistogramRegistry registry{"registry"};
  o2::analysis::rsn::Output* rsnOutput = nullptr;

  Service<o2::framework::O2DatabasePDG> pdg;

  int n = 0;
  float massPos = o2::track::PID::getMass(3);
  float massNeg = o2::track::PID::getMass(3);
  double* pointPair = nullptr;
  double* pointSys = nullptr;
  ROOT::Math::PxPyPzMVector d1, d2, mother;
  bool produceQA, dataQA, MCTruthQA, globalTrack, tpcPidOnly = false;
  float tpcnSigmaPos = 100.0f;
  float tpcnSigmaNeg = 100.0f;
  float combinedNSigma = 100.0f;
  float ptTOFThreshold = 0.5f;
  int tpcNClsFound = 70;
  int dauSize = 2;

  rsn::MixingType mixingType = rsn::MixingType::none;

  float vzCut = 10.0f;
  Filter triggerFilter = (o2::aod::evsel::sel8 == true);
  Filter vtxFilter = (nabs(o2::aod::collision::posZ) < vzCut);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>>;
  using EventCandidate = EventCandidates::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTOFFullKa>;

  using EventCandidatesMC = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>>;
  using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTOFFullKa, aod::McTrackLabels>;

  using EventCandidatesMCGen = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>;
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
    mixingType = rsn::mixingTypeName(static_cast<std::string>(produce.eventMixing));
    tpcnSigmaPos = static_cast<float>(cut.tpcnSigmaPos);
    tpcnSigmaNeg = static_cast<float>(cut.tpcnSigmaNeg);
    tpcNClsFound = static_cast<int>(cut.tpcNClsFound);
    globalTrack = static_cast<bool>(cut.globalTrack);
    combinedNSigma = static_cast<float>(cut.combinedNSigma);
    tpcPidOnly = static_cast<bool>(cut.tpcPidOnly);
    ptTOFThreshold = static_cast<float>(cut.ptTOFThreshold);

    pointPair = new double[static_cast<int>(o2::analysis::rsn::PairAxisType::unknown)];
    pointSys = new double[static_cast<int>(o2::analysis::rsn::SystematicsAxisType::unknown)];
    rsnOutput = new o2::analysis::rsn::OutputSparse();
    rsnOutput->init(sparseAxes, allAxes, sysAxes, allAxesSys, static_cast<bool>(produce.produceTrue), mixingType, static_cast<bool>(produce.produceLikesign), &registry);

    // Print summary of configuration
    LOGF(info, "=== PhianalysisTHnSparse configuration summary ===");
    LOGF(info, "produceQA: %s", produceQA ? "true" : "false");
    LOGF(info, "produceTrue: %s", static_cast<bool>(produce.produceTrue) ? "true" : "false");
    LOGF(info, "produceLikesign: %s", static_cast<bool>(produce.produceLikesign) ? "true" : "false");
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
    LOGF(info, "etapair: %.2f", static_cast<float>(cut.etapair));
    LOGF(info, "pt (min): %.2f", static_cast<float>(cut.pt));
    LOGF(info, "dcaXY: %.2f", static_cast<float>(cut.dcaXY));
    LOGF(info, "dcaZ: %.2f", static_cast<float>(cut.dcaZ));
    LOGF(info, "globalTrack: %s", globalTrack ? "true" : "false");
    LOGF(info, "tpcNClsFound: %d", tpcNClsFound);
    LOGF(info, "mixingType: %d", static_cast<int>(mixingType));
    LOGF(info, "numberofMixedEvents: %d", static_cast<int>(numberofMixedEvents));
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
      registry.add("QAEvent/hSelection", "Event selection statistics", kTH1D, {{3, 0.0f, 3.0f}});
      auto hEvent = registry.get<TH1>(HIST("QAEvent/hSelection"));
      hEvent->GetXaxis()->SetBinLabel(1, "Full event statistics");
      hEvent->GetXaxis()->SetBinLabel(2, "Events with at least one #phi candidate");
      hEvent->GetXaxis()->SetBinLabel(3, "#phi candidates");
      hEvent->SetMinimum(0.1);

      registry.add("QAEvent/hVtxZ", "Vertex position along the z-axis", kTH1F, {posZaxis});
      registry.add("QAEvent/hCent", "Distribution of multiplicity percentile", kTH1F, {{101, 0., 101.}});
      registry.add("QAEvent/hMult", "Multiplicity (amplitude of non-zero channels in the FV0A + FV0C) ", kTH1F, {{300, 0., 30000.}});

      // Track QA
      registry.add("QATrack/hSelection", "Tracks statistics", kTH1D, {{10, 0.0f, 10.0f}});
      auto hTrack = registry.get<TH1>(HIST("QATrack/hSelection"));
      hTrack->GetXaxis()->SetBinLabel(1, "all tracks");
      hTrack->GetXaxis()->SetBinLabel(2, "passed pT cut");
      hTrack->GetXaxis()->SetBinLabel(3, "passed eta cut");
      hTrack->GetXaxis()->SetBinLabel(4, "passed DCA cut");
      hTrack->GetXaxis()->SetBinLabel(5, "passed PID cut");
      hTrack->GetXaxis()->SetBinLabel(6, "passed rapidity cut");
      hTrack->GetXaxis()->SetBinLabel(7, "passed tpcNClsFound cut");
      hTrack->GetXaxis()->SetBinLabel(8, "passed isPrimaryTrack cut");
      hTrack->GetXaxis()->SetBinLabel(9, "passed isPVContributor cut");
      hTrack->GetXaxis()->SetBinLabel(10, "passed all cuts");
      hTrack->SetMinimum(0.1);

      registry.add("QATrack/hRapidity", "Rapidity distribution of K^{+} and K^{-}", kTH1F, {{200, -1, 1}});
      registry.add("QATrack/hEta", "Pseudorapidity distribution of K^{+} and K^{-}", kTH1F, {{200, -1, 1}});
      registry.add("QATrack/hTPCNClsFound", "Distribution of TPC NClsFound of K^{+} and K^{-}", kTH1F, {tpcNClsFoundQAaxis});
      registry.add("QATrack/hDCAxy", "Distribution of DCA_{xy} of K^{+} and K^{-}", kTH1F, {dcaXYaxis});
      registry.add("QATrack/hDCAz", "Distribution of DCA_{z} of K^{+} and K^{-}", kTH1F, {dcaZaxis});
      registry.add("QATrack/hPt", "Distribution of p_{T} of K^{+} and K^{-}", kTH1F, {ptaxis});

      // TPC PID
      registry.add("QATrack/hTPCnSigma", "Distribution of TPC nSigma of K^{+} and K^{-}", kTH1F, {{200, -10, 10}});

      registry.add("QATrack/h2TPCnSigma", "", kTH2F, {{200, -10, 10}, {200, -10, 10}});
      auto h2TPCnSigma = registry.get<TH2>(HIST("QATrack/h2TPCnSigma"));
      h2TPCnSigma->GetXaxis()->SetTitle("n#sigma_{TPC} K^{+}");
      h2TPCnSigma->GetYaxis()->SetTitle("n#sigma_{TPC} K^{-}");

      registry.add("QATrack/h2TPCnSigmaPt", "", kTH2F, {ptaxis, {200, -10, 10}});
      auto h2TPCnSigmaPt = registry.get<TH2>(HIST("QATrack/h2TPCnSigmaPt"));
      h2TPCnSigmaPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h2TPCnSigmaPt->GetYaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");

      // TOF PID
      registry.add("QATrack/hTOFnSigma", "Distribution of TOF nSigma of K^{+} and K^{-}", kTH1F, {{200, -10, 10}});

      registry.add("QATrack/h2TOFnSigma", "", kTH2F, {{200, -10, 10}, {200, -10, 10}});
      auto h2TOFnSigma = registry.get<TH2>(HIST("QATrack/h2TOFnSigma"));
      h2TOFnSigma->GetXaxis()->SetTitle("n#sigma_{TOF} K^{+}");
      h2TOFnSigma->GetYaxis()->SetTitle("n#sigma_{TOF} K^{-}");

      registry.add("QATrack/h2TOFnSigmaPt", "", kTH2F, {ptaxis, {200, -10, 10}});
      auto h2TOFnSigmaPt = registry.get<TH2>(HIST("QATrack/h2TOFnSigmaPt"));
      h2TOFnSigmaPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h2TOFnSigmaPt->GetYaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");

      // MC Truth
      if (static_cast<bool>(produce.produceTrue)) {
        registry.add("QAMC/Truth/hMCEvent", "MC Truth Event statistics", kTH1F, {{1, 0.0f, 1.0f}});
        auto hMCEventTruth = registry.get<TH1>(HIST("QAMC/Truth/hMCEvent"));
        hMCEventTruth->GetXaxis()->SetBinLabel(1, "Full MC Truth event statistics");
        hMCEventTruth->SetMinimum(0.1);

        registry.add("QAMC/hInvMassTrueFalse", "", kTH1F, {invAxis}); // not written events in True distribution due to repetition of mothers??

        // MC Gen
        registry.add("QAMC/Gen/hMCEvent", "MC Gen Event statistics", kTH1F, {{3, 0.0f, 3.0f}});
        auto hMCEventGen = registry.get<TH1>(HIST("QAMC/Gen/hMCEvent"));
        hMCEventGen->GetXaxis()->SetBinLabel(1, "Full McCollision statistics");
        hMCEventGen->GetXaxis()->SetBinLabel(2, "McCollision V_{z} cut");
        hMCEventGen->GetXaxis()->SetBinLabel(3, "collisions");
        hMCEventGen->SetMinimum(0.1);
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
    }

    pointSys[static_cast<int>(o2::analysis::rsn::SystematicsAxisType::ncl)] = tpcNClsFound;
    rsnOutput->fillSystematics(pointSys);
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
      if (std::abs(track.tpcNSigmaKa()) >= nSigmaCut)
        return false;
    } else {
      if (std::sqrt(track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa()) >= combinedNSigma)
        return false;
    }
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 4.5);

    // Apply rapidity cut
    if (std::abs(track.rapidity(isPositive ? massPos : massNeg)) > static_cast<float>(cut.rapidity)) {
      return false;
    }
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 5.5);

    // Apply tpcNClsFound cut
    if (track.tpcNClsFound() < tpcNClsFound)
      return false;
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 6.5);

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
      registry.fill(HIST("QATrack/hSelection"), 7.5);

    // Apply PV Contributor cuts
    if (!track.isPVContributor())
      return false;
    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 8.5);

    if (produceQA && dataQA)
      registry.fill(HIST("QATrack/hSelection"), 9.5);

    return true;
  }
  template <typename T>
  bool selectedPair(ROOT::Math::PxPyPzMVector& mother, const T& track1, const T& track2)
  {
    d1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPos);
    d2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massNeg);
    mother = d1 + d2;

    if (std::abs(mother.Eta()) > static_cast<float>(cut.etapair))
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
    n = 0;
    if (produceQA) {
      registry.fill(HIST("QAEvent/hSelection"), 0.5);
      registry.fill(HIST("QAEvent/hVtxZ"), collision.posZ());
      registry.fill(HIST("QAEvent/hMult"), getMultiplicity(collision));
      registry.fill(HIST("QAEvent/hCent"), getCentrality(collision));

      dataQA = true;
      for (const auto& track : posDaughters) {
        selectedTrack(track, true);
      }
      for (const auto& track : negDaughters) {
        selectedTrack(track, false);
      }
      dataQA = false;
    }

    if (static_cast<int>(verbose.verboselevel) > 0 && static_cast<int>(verbose.refresh) > 0 && collision.globalIndex() % static_cast<int>(verbose.refresh) == static_cast<int>(verbose.refreshIndex))
      LOGF(info, "%d pos=%lld neg=%lld, Z vertex position: %f [cm]", collision.globalIndex(), posDaughters.size(), negDaughters.size(), collision.posZ());

    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughters, negDaughters))) {

      if (!selectedTrack(track1, true)) // track1 is positive
        continue;
      if (!selectedTrack(track2, false)) // track2 is negative
        continue;

      if (!selectedPair(mother, track1, track2))
        continue;

      if (produceQA) {
        registry.fill(HIST("QATrack/h2TPCnSigma"), track1.tpcNSigmaKa(), track2.tpcNSigmaKa());
        registry.fill(HIST("QATrack/h2TPCnSigmaPt"), track1.pt(), track1.tpcNSigmaKa());
        registry.fill(HIST("QATrack/h2TPCnSigmaPt"), track2.pt(), track2.tpcNSigmaKa());

        registry.fill(HIST("QATrack/h2TOFnSigma"), track1.tofNSigmaKa(), track2.tofNSigmaKa());
        registry.fill(HIST("QATrack/h2TOFnSigmaPt"), track1.pt(), track1.tofNSigmaKa());
        registry.fill(HIST("QATrack/h2TOFnSigmaPt"), track2.pt(), track2.tofNSigmaKa());

        registry.fill(HIST("QATrack/hTPCnSigma"), track1.tpcNSigmaKa());
        registry.fill(HIST("QATrack/hTPCnSigma"), track2.tpcNSigmaKa());
        if (track1.hasTOF())
          registry.fill(HIST("QATrack/hTOFnSigma"), track1.tofNSigmaKa());
        if (track2.hasTOF())
          registry.fill(HIST("QATrack/hTOFnSigma"), track2.tofNSigmaKa());

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
      }

      pointPair = fillPointPair(mother.M(),
                                mother.Pt(),
                                getMultiplicity(collision),
                                getCentrality(collision),
                                track1.tpcNSigmaKa(),
                                track2.tpcNSigmaKa(),
                                mother.Eta(),
                                mother.Rapidity(),
                                collision.posZ(),
                                0,
                                0,
                                0);
      rsnOutput->fillUnlikepm(pointPair);

      if (produceQA)
        registry.fill(HIST("QAEvent/hSelection"), 2.5);
      if (n == 0)
        registry.fill(HIST("QAEvent/hSelection"), 1.5);
      n = n + 1;
    }

    if (static_cast<bool>(produce.produceLikesign)) {

      for (const auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posDaughters, posDaughters))) {
        if (!selectedTrack(track1, true)) // both positive
          continue;
        if (!selectedTrack(track2, true)) // both positive
          continue;

        if (!selectedPair(mother, track1, track2))
          continue;

        if (static_cast<int>(verbose.verboselevel) > 1)
          LOGF(info, "Like-sign positive: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.M());

        pointPair = fillPointPair(mother.M(),
                                  mother.Pt(),
                                  getMultiplicity(collision),
                                  getCentrality(collision),
                                  track1.tpcNSigmaKa(),
                                  track2.tpcNSigmaKa(),
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

        if (!selectedPair(mother, track1, track2))
          continue;

        if (static_cast<int>(verbose.verboselevel) > 1)
          LOGF(info, "Like-sign negative: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.M());

        pointPair = fillPointPair(mother.M(),
                                  mother.Pt(),
                                  getMultiplicity(collision),
                                  getCentrality(collision),
                                  track1.tpcNSigmaKa(),
                                  track2.tpcNSigmaKa(),
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

  void processTrue(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& /*tracks*/, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    if (!static_cast<bool>(produce.produceTrue))
      return;

    registry.fill(HIST("QAMC/Truth/hMCEvent"), 0.5);

    auto posDaughtersMC = positiveMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDaughtersMC = negativeMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (!collision.has_mcCollision()) {
      if (static_cast<int>(verbose.verboselevel) > 0)
        LOGF(warning, "No MC collision for this collision, skip...");
      return;
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

          if (static_cast<int>(verbose.verboselevel) > 1) {
            LOGF(info, "True: %d, d1=%d (%ld), d2=%d (%ld), mother=%d (%ld)", n, mctrack1.pdgCode(), mctrack1.globalIndex(), mctrack2.pdgCode(), mctrack2.globalIndex(), mothertrack1.pdgCode(), mothertrack1.globalIndex());
            LOGF(info, "%d px: %f, py=%f, pz=%f, px: %f, py=%f, pz=%f", n, mctrack1.px(), mctrack1.py(), mctrack1.pz(), mctrack2.px(), mctrack2.py(), mctrack2.pz());
          }

          if (!selectedPair(mother, mctrack1, mctrack2))
            continue;

          if (n > 0) {
            if (produceQA)
              registry.fill(HIST("QAMC/hInvMassTrueFalse"), mother.M());
            continue;
          }

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    getMultiplicity(collision),
                                    getCentrality(collision),
                                    track1.tpcNSigmaKa(),
                                    track2.tpcNSigmaKa(),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    collision.posZ(),
                                    0,
                                    0,
                                    0);

          rsnOutput->fillUnliketrue(pointPair);
          n++;
        }
      }
    }
  }
  PROCESS_SWITCH(PhianalysisTHnSparse, processTrue, "Process Event for MC reconstruction.", false);

  void processGen(aod::McCollision const& mcCollision, soa::SmallGroups<EventCandidatesMCGen> const& collisions, LabeledTracks const& /*particles*/, aod::McParticles const& mcParticles)
  {

    if (produceQA)
      registry.fill(HIST("QAMC/Gen/hMCEvent"), 0.5);

    if (std::abs(mcCollision.posZ()) > vzCut)
      return;

    if (produceQA)
      registry.fill(HIST("QAMC/Gen/hMCEvent"), 1.5);

    if (collisions.size() == 0)
      return;

    for (const auto& collision : collisions) {
      if (produceQA)
        registry.fill(HIST("QAMC/Gen/hMCEvent"), 2.5);

      if (!collision.has_mcCollision()) {
        if (static_cast<int>(verbose.verboselevel) > 0)
          LOGF(warning, "No McCollision for this collision, skip...");
        return;
      }

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
        }
      }
    }
  }
  PROCESS_SWITCH(PhianalysisTHnSparse, processGen, "Process MC Generated.", false);

  void processMixed(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    if (mixingType == rsn::MixingType::none)
      return;

    auto tracksTuple = std::make_tuple(tracks);

    BinningTypeVzCe binningVzCe{{axisVertexMixing, axisCentralityMixing}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVzCe> pairVzCe{binningVzCe, numberofMixedEvents, -1, collisions, tracksTuple, &cache};

    BinningTypeVzMu binningVzMu{{axisVertexMixing, axisMultiplicityMixing}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVzMu> pairVzMu{binningVzMu, numberofMixedEvents, -1, collisions, tracksTuple, &cache};

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

          if (!selectedPair(mother, track1, track2))
            continue;

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    track1.tpcNSigmaKa(),
                                    track2.tpcNSigmaKa(),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    c1.posZ(),
                                    getMultiplicity(c2),
                                    getCentrality(c2),
                                    c2.posZ());

          rsnOutput->fillMixingpm(pointPair);
        }

        if (static_cast<bool>(produce.produceLikesign)) {

          for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersc1, posDaughtersc2))) {

            if (!selectedTrack(track1, true)) // track1 is positive
              continue;
            if (!selectedTrack(track2, true)) // track2 is positive
              continue;

            if (!selectedPair(mother, track1, track2))
              continue;

            pointPair = fillPointPair(mother.M(),
                                      mother.Pt(),
                                      getMultiplicity(c1),
                                      getCentrality(c1),
                                      track1.tpcNSigmaKa(),
                                      track2.tpcNSigmaKa(),
                                      mother.Eta(),
                                      mother.Rapidity(),
                                      c1.posZ(),
                                      getMultiplicity(c2),
                                      getCentrality(c2),
                                      c2.posZ());

            rsnOutput->fillMixingpp(pointPair);
          }

          for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(negDaughtersc1, negDaughtersc2))) {

            if (!selectedTrack(track1, false))
              continue;
            if (!selectedTrack(track2, false))
              continue;

            if (!selectedPair(mother, track1, track2))
              continue;
            pointPair = fillPointPair(mother.M(),
                                      mother.Pt(),
                                      getMultiplicity(c1),
                                      getCentrality(c1),
                                      track1.tpcNSigmaKa(),
                                      track2.tpcNSigmaKa(),
                                      mother.Eta(),
                                      mother.Rapidity(),
                                      c1.posZ(),
                                      getMultiplicity(c2),
                                      getCentrality(c2),
                                      c2.posZ());

            rsnOutput->fillMixingmm(pointPair);
          }
        }

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersc2, negDaughtersc1))) {

          if (!selectedTrack(track1, true)) // track1 is positive
            continue;
          if (!selectedTrack(track2, false)) // track2 is negative
            continue;

          if (!selectedPair(mother, track1, track2))
            continue;

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    track1.tpcNSigmaKa(),
                                    track2.tpcNSigmaKa(),
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

          if (!selectedPair(mother, track1, track2))
            continue;

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    track1.tpcNSigmaKa(),
                                    track2.tpcNSigmaKa(),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    c1.posZ(),
                                    getMultiplicity(c2),
                                    getCentrality(c2),
                                    c2.posZ());

          rsnOutput->fillMixingpm(pointPair);
        }

        if (static_cast<bool>(produce.produceLikesign)) {

          for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersc1, posDaughtersc2))) {

            if (!selectedTrack(track1, true)) // track1 is positive

              continue;
            if (!selectedTrack(track2, true)) // track2 is positive
              continue;

            if (!selectedPair(mother, track1, track2))
              continue;

            pointPair = fillPointPair(mother.M(),
                                      mother.Pt(),
                                      getMultiplicity(c1),
                                      getCentrality(c1),
                                      track1.tpcNSigmaKa(),
                                      track2.tpcNSigmaKa(),
                                      mother.Eta(),
                                      mother.Rapidity(),
                                      c1.posZ(),
                                      getMultiplicity(c2),
                                      getCentrality(c2),
                                      c2.posZ());

            rsnOutput->fillMixingpp(pointPair);
          }

          for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(negDaughtersc1, negDaughtersc2))) {

            if (!selectedTrack(track1, false))

              continue;
            if (!selectedTrack(track2, false))
              continue;

            if (!selectedPair(mother, track1, track2))
              continue;
            pointPair = fillPointPair(mother.M(),
                                      mother.Pt(),
                                      getMultiplicity(c1),
                                      getCentrality(c1),
                                      track1.tpcNSigmaKa(),
                                      track2.tpcNSigmaKa(),
                                      mother.Eta(),
                                      mother.Rapidity(),
                                      c1.posZ(),
                                      getMultiplicity(c2),
                                      getCentrality(c2),
                                      c2.posZ());

            rsnOutput->fillMixingmm(pointPair);
          }
        }

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDaughtersc2, negDaughtersc1))) {

          if (!selectedTrack(track1, true))

            continue;
          if (!selectedTrack(track2, false))
            continue;

          if (!selectedPair(mother, track1, track2))
            continue;

          pointPair = fillPointPair(mother.M(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    track1.tpcNSigmaKa(),
                                    track2.tpcNSigmaKa(),
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
