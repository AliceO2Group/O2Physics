// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
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

#include <TLorentzVector.h>
#include <vector>
#include <string>

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "PWGLF/Utils/rsnOutput.h"
// #include "TDatabasePDG.h"

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

  struct : ConfigurableGroup {
    Configurable<int> verboselevel{"verboselevel", 0, "Verbose level"};
    Configurable<int> refresh{"refresh", 0, "Freqency of print event information."};
    Configurable<int> refreshIndex{"refreshIndex", 0, "Freqency of print event information index."};
  } verbose;

  Configurable<int> dautherPos{"dautherPos", 3, "Particle type of the positive dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> dautherNeg{"dautherNeg", 3, "Particle type of the negative dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> motherPDG{"motherPDG", 333, "PDG code of mother particle."};
  Configurable<int> dautherPosPDG{"dautherPosPDG", 321, "PDG code of positive dauther particle."};
  Configurable<int> dautherNegPDG{"dautherNegPDG", 321, "PDG code of negative dauther particle."};

  struct : ConfigurableGroup {
    Configurable<float> tpcnSigmaPos{"tpcnSigmaPos", 3.0f, "TPC NSigma cut of the positive particle."};
    Configurable<float> tpcnSigmaNeg{"tpcnSigmaNeg", 3.0f, "TPC NSigma cut of the negative particle."};
    Configurable<float> vZ{"vZ", 10.0f, "Z vertex range."};
    Configurable<float> y{"y", 0.5, "Rapidity cut (maximum)."};
    Configurable<float> etatrack{"etatrack", 0.8, "Eta cut for track."};
    Configurable<float> etapair{"etapair", 0.8, "Eta cut for pair."};
    Configurable<float> pt{"pt", 0.15f, "Cut: Minimal value of tracks pt."};
    Configurable<float> dcaXY{"dcaXY", 1.0f, "Cut: Maximal value of tracks DCA XY."};
    Configurable<float> dcaZ{"dcaZ", 1.0f, "Cut: Maximal value of tracks DCA Z."};
    Configurable<int> tpcNClsFound{"tpcNClsFound", 70, "Cut: Minimal value of found TPC clasters"};
  } cut;

  Configurable<std::vector<std::string>> sparseAxes{"sparseAxes", std::vector<std::string>{o2::analysis::rsn::pair_axis::names}, "Axes."};
  Configurable<std::vector<std::string>> sysAxes{"sysAxes", std::vector<std::string>{o2::analysis::rsn::systematic_axis::names}, "Axes."};

  ConfigurableAxis invaxis{"invaxis", {130, 0.97, 1.1}, "Invariant mass axis binning."};
  ConfigurableAxis ptaxis{"ptaxis", {20, 0., 20.}, "Pt axis binning."};
  ConfigurableAxis vzaxis{"vzaxis", {40, -20., 20.}, "Z vertex position axis binning."};
  ConfigurableAxis multiplicityaxis{"multiplicityaxis", {50, 0., 5000.}, "Multiplicity axis binning."};
  ConfigurableAxis centralityaxis{"centralityaxis", {20, 0., 100.}, "Centrality axis binning."};
  ConfigurableAxis etaaxis{"etaaxis", {16., -1.0 * static_cast<float>(cut.etatrack), static_cast<float>(cut.etatrack)}, "Pseudorapidity axis binning."};
  ConfigurableAxis rapidityaxis{"rapidityaxis", {10., -1.0 * static_cast<float>(cut.y), static_cast<float>(cut.y)}, "Rapidity axis binning."};
  ConfigurableAxis nsigmaaxisPos{"nsigmaaxisPos", {1, 0., static_cast<float>(cut.tpcnSigmaPos)}, "NSigma of positive particle axis binning in THnSparse."};
  ConfigurableAxis nsigmaaxisNeg{"nsigmaaxisNeg", {1, 0., static_cast<float>(cut.tpcnSigmaNeg)}, "NSigma of negative particle axis binning in THnSparse."};

  // mixing
  using BinningTypeVzMu = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
  using BinningTypeVzCe = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  Configurable<int> numberofMixedEvents{"numberofMixedEvents", 5, "Number of events that should be mixed."};
  ConfigurableAxis axisVertexMixing{"axisVertexMixing", {5, -10, 10}, "Z vertex axis binning for mixing"};
  ConfigurableAxis axisMultiplicityMixing{"axisMultiplicityMixing", {5, 0, 5000}, "TPC multiplicity for bin"};
  ConfigurableAxis axisCentralityMixing{"axisCentralityMixing", {10, 0, 100}, "Multiplicity percentil binning for mixing"};

  // defined in DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h
  float massPos = o2::track::PID::getMass(dautherPos);
  float massNeg = o2::track::PID::getMass(dautherNeg);

  // Axes specifications
  AxisSpec posZaxis = {400, -20., 20., "V_{z} (cm)"};
  AxisSpec dcaXYaxis = {120, -3.0, 3.0, "DCA_{xy} (cm)"};
  AxisSpec dcaZaxis = {120, -3.0, 3.0, "DCA_{z} (cm)"};
  AxisSpec etaQAaxis = {1000, -1.0, 1.0, "#eta"};
  AxisSpec tpcNClsFoundQAaxis = {110, 50., 160., "tpcNClsFound"};

  HistogramRegistry registry{"registry"};
  o2::analysis::rsn::Output* rsnOutput = nullptr;

  Service<o2::framework::O2DatabasePDG> pdg;

  int n = 0;
  double* pointPair = nullptr;
  double* pointSys = nullptr;
  TLorentzVector d1, d2, mother;
  bool produceQA, dataQA, MCTruthQA, t1, t2 = false;
  int id;
  float etapair = 0.8;
  float tpcnSigmaPos = 3.0f;
  float tpcnSigmaNeg = 3.0f;
  int tpcNClsFound = 70;

  rsn::MixingType mixingType = rsn::MixingType::none;

  Filter triggerFilter = (o2::aod::evsel::sel8 == true);
  Filter vtxFilter = (nabs(o2::aod::collision::posZ) < static_cast<float>(cut.vZ));

  Filter ptFilter = nabs(aod::track::pt) > static_cast<float>(cut.pt);
  Filter etaFilter = nabs(aod::track::eta) < static_cast<float>(cut.etatrack);
  Filter dcaFilter = (nabs(o2::aod::track::dcaXY) < static_cast<float>(cut.dcaXY)) && (nabs(o2::aod::track::dcaZ) < static_cast<float>(cut.dcaZ));

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>>;
  using EventCandidate = EventCandidates::iterator;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa>>;

  using EventCandidatesMC = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::McTrackLabels>>;

  using EventCandidatesMCGen = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>;

  using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  Partition<TrackCandidates> positive = (aod::track::signed1Pt > 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < std::abs(static_cast<float>(cut.tpcnSigmaPos)));
  Partition<TrackCandidates> negative = (aod::track::signed1Pt < 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < std::abs(static_cast<float>(cut.tpcnSigmaNeg)));

  Partition<TrackCandidatesMC> positiveMC = (aod::track::signed1Pt > 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < std::abs(static_cast<float>(cut.tpcnSigmaPos)));
  Partition<TrackCandidatesMC> negativeMC = (aod::track::signed1Pt < 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < std::abs(static_cast<float>(cut.tpcnSigmaNeg)));

  void init(o2::framework::InitContext&)
  {
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
    etapair = static_cast<float>(cut.etapair);
    tpcnSigmaPos = static_cast<float>(cut.tpcnSigmaPos);
    tpcnSigmaNeg = static_cast<float>(cut.tpcnSigmaNeg);
    tpcNClsFound = static_cast<int>(cut.tpcNClsFound);

    pointPair = new double[static_cast<int>(o2::analysis::rsn::PairAxisType::unknown)];
    pointSys = new double[static_cast<int>(o2::analysis::rsn::SystematicsAxisType::unknown)];
    rsnOutput = new o2::analysis::rsn::OutputSparse();
    rsnOutput->init(sparseAxes, allAxes, sysAxes, allAxesSys, static_cast<bool>(produce.produceTrue), mixingType, static_cast<bool>(produce.produceLikesign), &registry);

    if (produceQA) {
      // Event QA
      registry.add("QAEvent/hSelection", "Event selection statistics", kTH1F, {{1, 0.0f, 1.0f}});
      auto hEvent = registry.get<TH1>(HIST("QAEvent/hSelection"));
      hEvent->GetXaxis()->SetBinLabel(1, "Full event statistics");
      hEvent->SetMinimum(0.1);

      registry.add("QAEvent/hVtxZ", "Vertex position along the z-axis", kTH1F, {posZaxis});
      registry.add("QAEvent/hCent", "Distribution of multiplicity percentile", kTH1F, {{100, 0., 100.}});
      registry.add("QAEvent/hMult", "Multiplicity (amplitude of non-zero channels in the FV0A + FV0C) ", kTH1F, {{300, 0., 30000.}});
      registry.add("QAEvent/h2Size", "Number of positive vs. negative Kaons per collision", kTH2F, {{30, 0., 30.}, {30, 0., 30.}});

      // Track QA
      registry.add("QATrack/hSelection", "Tracks combinations statistics", kTH1F, {{11, 0.0f, 11.0f}});
      auto hTrack = registry.get<TH1>(HIST("QATrack/hSelection"));
      hTrack->GetXaxis()->SetBinLabel(1, "all K^{+} K^{-} combinations");
      hTrack->GetXaxis()->SetBinLabel(2, "all K^{+}");
      hTrack->GetXaxis()->SetBinLabel(3, "all K^{-}");
      hTrack->GetXaxis()->SetBinLabel(4, "K^{+} tpcNClsFound");
      hTrack->GetXaxis()->SetBinLabel(5, "K^{-} tpcNClsFound");
      hTrack->GetXaxis()->SetBinLabel(6, "K^{+} isPrimaryTrack");
      hTrack->GetXaxis()->SetBinLabel(7, "K^{-} isPrimaryTrack");
      hTrack->GetXaxis()->SetBinLabel(8, "K^{+} isPVContributor");
      hTrack->GetXaxis()->SetBinLabel(9, "K^{-} isPVContributor");
      hTrack->GetXaxis()->SetBinLabel(10, "selected combinations");
      hTrack->GetXaxis()->SetBinLabel(11, "selected pairs (eta cut)");
      hTrack->SetMinimum(0.1);

      // Track1
      registry.add("QATrack/bs/hTrack1pt", "K^{+} p_{T} before selection", kTH1F, {ptAxis});
      registry.add("QATrack/bs/hTrack1DCAxy", "K^{+} DCA_{xy} before selection", kTH1F, {dcaXYaxis});
      registry.add("QATrack/bs/hTrack1DCAz", "K^{+} DCA_{z} before selection", kTH1F, {dcaZaxis});
      registry.add("QATrack/bs/hTrack1eta", "K^{+} #eta before selection", kTH1F, {{etaQAaxis}});
      registry.add("QATrack/bs/hTrack1tpcNClsFound", "K^{+} tpcNClsFound before selection", kTH1F, {{tpcNClsFoundQAaxis}});

      registry.add("QATrack/as/hTrack1pt", "K^{+} p_{T} after selection", kTH1F, {ptAxis});
      registry.add("QATrack/as/hTrack1DCAxy", "K^{+} DCA_{xy} after selection", kTH1F, {dcaXYaxis});
      registry.add("QATrack/as/hTrack1DCAz", "K^{+} DCA_{z} after selection", kTH1F, {dcaZaxis});
      registry.add("QATrack/as/hTrack1eta", "K^{+} #eta after selection", kTH1F, {{etaQAaxis}});
      registry.add("QATrack/as/hTrack1tpcNClsFound", "K^{+} tpcNClsFound after selection", kTH1F, {{tpcNClsFoundQAaxis}});

      // Track2
      registry.add("QATrack/bs/hTrack2pt", "K^{-} p_{T} before selection", kTH1F, {ptAxis});
      registry.add("QATrack/bs/hTrack2DCAxy", "K^{-} DCA_{xy} before selection", kTH1F, {dcaXYaxis});
      registry.add("QATrack/bs/hTrack2DCAz", "K^{-} DCA_{z} before selection", kTH1F, {dcaZaxis});
      registry.add("QATrack/bs/hTrack2eta", "K^{-} #eta before selection", kTH1F, {{etaQAaxis}});
      registry.add("QATrack/bs/hTrack2tpcNClsFound", "K^{-} tpcNClsFound before selection", kTH1F, {{tpcNClsFoundQAaxis}});

      registry.add("QATrack/as/hTrack2pt", "K^{-} p_{T} after selection", kTH1F, {ptAxis});
      registry.add("QATrack/as/hTrack2DCAxy", "K^{-} DCA_{xy} after selection", kTH1F, {dcaXYaxis});
      registry.add("QATrack/as/hTrack2DCAz", "K^{-} DCA_{z} after selection", kTH1F, {dcaZaxis});
      registry.add("QATrack/as/hTrack2eta", "K^{-} #eta after selection", kTH1F, {{etaQAaxis}});
      registry.add("QATrack/as/hTrack2tpcNClsFound", "K^{-} tpcNClsFound after selection", kTH1F, {{tpcNClsFoundQAaxis}});

      // TPC PID
      registry.add("QATrack/TPCPID/h2TracknSigma", "", kTH2F, {{120, -6, 6}, {120, -6, 6}});
      auto h2TracknSigma = registry.get<TH2>(HIST("QATrack/TPCPID/h2TracknSigma"));
      h2TracknSigma->GetXaxis()->SetTitle("n#sigma_{TPC} K^{+}");
      h2TracknSigma->GetYaxis()->SetTitle("n#sigma_{TPC} K^{-}");

      registry.add("QATrack/TPCPID/h2nTrack1SigmaPt", "", kTH2F, {{100, 0, 10}, {120, -6, 6}});
      auto h2nTrack1SigmaPt = registry.get<TH2>(HIST("QATrack/TPCPID/h2nTrack1SigmaPt"));
      h2nTrack1SigmaPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h2nTrack1SigmaPt->GetYaxis()->SetTitle("n#sigma_{TPC} K^{+}");

      registry.add("QATrack/TPCPID/h2nTrack2SigmaPt", "", kTH2F, {{100, 0, 10}, {120, -6, 6}});
      auto h2nTrack2SigmaPt = registry.get<TH2>(HIST("QATrack/TPCPID/h2nTrack2SigmaPt"));
      h2nTrack2SigmaPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h2nTrack2SigmaPt->GetYaxis()->SetTitle("n#sigma_{TPC} K^{-}");

      // MC Truth
      if (static_cast<bool>(produce.produceTrue)) {
        registry.add("QAMC/Truth/hMCEvent", "MC Truth Event statistics", kTH1F, {{1, 0.0f, 1.0f}});
        auto hMCEventTruth = registry.get<TH1>(HIST("QAMC/Truth/hMCEvent"));
        hMCEventTruth->GetXaxis()->SetBinLabel(1, "Full MC Truth event statistics");
        hMCEventTruth->SetMinimum(0.1);

        registry.add("QAMC/Truth/hMCTrack", "MC Truth Track statistics", kTH1F, {{17, 0.0f, 17.0f}});
        auto hMCTrackTruth = registry.get<TH1>(HIST("QAMC/Truth/hMCTrack"));
        hMCTrackTruth->GetXaxis()->SetBinLabel(1, "all K^{+} K^{-} combinations");
        hMCTrackTruth->GetXaxis()->SetBinLabel(2, "all K^{+}");
        hMCTrackTruth->GetXaxis()->SetBinLabel(3, "all K^{-}");
        hMCTrackTruth->GetXaxis()->SetBinLabel(4, "K^{+} tpcNClsFound");
        hMCTrackTruth->GetXaxis()->SetBinLabel(5, "K^{-} tpcNClsFound");
        hMCTrackTruth->GetXaxis()->SetBinLabel(6, "K^{+} isPrimaryTrack");
        hMCTrackTruth->GetXaxis()->SetBinLabel(7, "K^{-} isPrimaryTrack");
        hMCTrackTruth->GetXaxis()->SetBinLabel(8, "K^{+} isPVContributor");
        hMCTrackTruth->GetXaxis()->SetBinLabel(9, "K^{-} isPVContributor");
        hMCTrackTruth->GetXaxis()->SetBinLabel(10, "selected combinations");
        hMCTrackTruth->GetXaxis()->SetBinLabel(11, "MCtrack PDG = 321");
        hMCTrackTruth->GetXaxis()->SetBinLabel(12, "all mothers");
        hMCTrackTruth->GetXaxis()->SetBinLabel(13, "equal mother PDGs");
        hMCTrackTruth->GetXaxis()->SetBinLabel(14, "equal mother IDs");
        hMCTrackTruth->GetXaxis()->SetBinLabel(15, "mother rapidity cut");
        hMCTrackTruth->GetXaxis()->SetBinLabel(16, "mother PDG = 333");
        hMCTrackTruth->GetXaxis()->SetBinLabel(17, "selected pairs (eta cut)");
        hMCTrackTruth->SetMinimum(0.1);

        registry.add("QAMC/hInvMassTrueFalse", "", kTH1F, {invAxis}); // not written events in True distribution due to repetition of mothers??

        // MC Gen
        registry.add("QAMC/Gen/hMCEvent", "MC Gen Event statistics", kTH1F, {{3, 0.0f, 3.0f}});
        auto hMCEventGen = registry.get<TH1>(HIST("QAMC/Gen/hMCEvent"));
        hMCEventGen->GetXaxis()->SetBinLabel(1, "Full McCollision statistics");
        hMCEventGen->GetXaxis()->SetBinLabel(2, "McCollision V_{z} cut");
        hMCEventGen->GetXaxis()->SetBinLabel(3, "collisions");
        hMCEventGen->SetMinimum(0.1);

        registry.add("QAMC/Gen/hMCTrack", "MC Gen Track statistics", kTH1D, {{7, 0.0f, 7.0f}});
        auto hMCTrackGen = registry.get<TH1>(HIST("QAMC/Gen/hMCTrack"));
        hMCTrackGen->GetXaxis()->SetBinLabel(1, "all mcParticles");
        hMCTrackGen->GetXaxis()->SetBinLabel(2, "rapidity cut");
        hMCTrackGen->GetXaxis()->SetBinLabel(3, "particle PDG = 333");
        hMCTrackGen->GetXaxis()->SetBinLabel(4, "has 2 dauthers");
        hMCTrackGen->GetXaxis()->SetBinLabel(5, "all dauthers");
        hMCTrackGen->GetXaxis()->SetBinLabel(6, "isPhysicalPrimary");
        hMCTrackGen->GetXaxis()->SetBinLabel(7, "selected pairs");
        hMCTrackGen->SetMinimum(0.1);
      }
      // Mixing QA
      if (mixingType != rsn::MixingType::none) {
        registry.add("QAMixing/hSelection", "Event mixing selection statistics", kTH1F, {{1, 0.0f, 1.0f}});
        auto hEM = registry.get<TH1>(HIST("QAMixing/hSelection"));
        hEM->GetXaxis()->SetBinLabel(1, "Full event mixing statistics");
        hEM->SetMinimum(0.1);

        registry.add("QAMixing/hTrackSelection", "Event mixing tracks combinations statistics", kTH1F, {{4, 0.0f, 4.0f}});
        auto hEMTrack = registry.get<TH1>(HIST("QAMixing/hTrackSelection"));
        hEMTrack->GetXaxis()->SetBinLabel(1, "all K^{+} K^{-} combinations");
        hEMTrack->GetXaxis()->SetBinLabel(2, "all K^{+}");
        hEMTrack->GetXaxis()->SetBinLabel(3, "all K^{-}");
        hEMTrack->GetXaxis()->SetBinLabel(4, "selected pairs (eta cut)");
        hEMTrack->SetMinimum(0.1);

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
  bool selectedTrack(const T& track)
  {
    if (produceQA) {
      if (t1) {
        if (dataQA)
          registry.fill(HIST("QATrack/hSelection"), 1.5);
        if (MCTruthQA)
          registry.fill(HIST("QAMC/Truth/hMCTrack"), 1.5);
      }
      if (t2) {
        if (dataQA)
          registry.fill(HIST("QATrack/hSelection"), 2.5);
        if (MCTruthQA)
          registry.fill(HIST("QAMC/Truth/hMCTrack"), 2.5);
      }
    }

    if (track.tpcNClsFound() < tpcNClsFound)
      return false;
    if (produceQA) {
      if (t1) {
        if (dataQA)
          registry.fill(HIST("QATrack/hSelection"), 3.5);
        if (MCTruthQA)
          registry.fill(HIST("QAMC/Truth/hMCTrack"), 3.5);
      }
      if (t2) {
        if (dataQA)
          registry.fill(HIST("QATrack/hSelection"), 4.5);
        if (MCTruthQA)
          registry.fill(HIST("QAMC/Truth/hMCTrack"), 4.5);
      }
    }
    if (!track.isPrimaryTrack())
      return false;
    if (produceQA) {
      if (t1) {
        if (dataQA)
          registry.fill(HIST("QATrack/hSelection"), 5.5);
        if (MCTruthQA)
          registry.fill(HIST("QAMC/Truth/hMCTrack"), 5.5);
      }
      if (t2) {
        if (dataQA)
          registry.fill(HIST("QATrack/hSelection"), 6.5);
        if (MCTruthQA)
          registry.fill(HIST("QAMC/Truth/hMCTrack"), 6.5);
      }
    }
    if (!track.isPVContributor())
      return false;
    if (produceQA) {
      if (t1) {
        if (dataQA)
          registry.fill(HIST("QATrack/hSelection"), 7.5);
        if (MCTruthQA)
          registry.fill(HIST("QAMC/Truth/hMCTrack"), 7.5);
      }
      if (t2) {
        if (dataQA)
          registry.fill(HIST("QATrack/hSelection"), 8.5);
        if (MCTruthQA)
          registry.fill(HIST("QAMC/Truth/hMCTrack"), 8.5);
      }
    }
    return true;
  }
  template <typename T>
  bool selectedPair(TLorentzVector& mother, const T& track1, const T& track2)
  {
    d1.SetXYZM(track1.px(), track1.py(), track1.pz(), massPos);
    d2.SetXYZM(track2.px(), track2.py(), track2.pz(), massNeg);
    mother = d1 + d2;

    if (std::abs(mother.Eta()) > etapair)
      return false;
    return true;
  }
  template <typename T>
  float getMultiplicity(const T& collision)
  {
    float multiplicity = collision.multFT0C() + collision.multFT0A();
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
    auto posDauthers = positive->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDauthers = negative->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (produceQA) {
      registry.fill(HIST("QAEvent/hSelection"), 0.5);
      registry.fill(HIST("QAEvent/h2Size"), posDauthers.size(), negDauthers.size());
      registry.fill(HIST("QAEvent/hVtxZ"), collision.posZ());
      registry.fill(HIST("QAEvent/hMult"), getMultiplicity(collision));
      registry.fill(HIST("QAEvent/hCent"), getCentrality(collision));
    }

    if (static_cast<int>(verbose.verboselevel) > 0 && static_cast<int>(verbose.refresh) > 0 && collision.globalIndex() % static_cast<int>(verbose.refresh) == static_cast<int>(verbose.refreshIndex))
      LOGF(info, "%d pos=%lld neg=%lld, Z vertex position: %f [cm]", collision.globalIndex(), posDauthers.size(), negDauthers.size(), collision.posZ());

    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthers, negDauthers))) {
      if (produceQA) {
        registry.fill(HIST("QATrack/hSelection"), 0.5);

        registry.fill(HIST("QATrack/bs/hTrack1pt"), track1.pt());
        registry.fill(HIST("QATrack/bs/hTrack1DCAxy"), track1.dcaXY());
        registry.fill(HIST("QATrack/bs/hTrack1DCAz"), track1.dcaZ());
        registry.fill(HIST("QATrack/bs/hTrack1eta"), track1.eta());
        registry.fill(HIST("QATrack/bs/hTrack1tpcNClsFound"), track1.tpcNClsFound());

        registry.fill(HIST("QATrack/bs/hTrack2pt"), track2.pt());
        registry.fill(HIST("QATrack/bs/hTrack2DCAxy"), track2.dcaXY());
        registry.fill(HIST("QATrack/bs/hTrack2DCAz"), track2.dcaZ());
        registry.fill(HIST("QATrack/bs/hTrack2eta"), track2.eta());
        registry.fill(HIST("QATrack/bs/hTrack2tpcNClsFound"), track2.tpcNClsFound());
      }

      dataQA = true;
      t1 = true;
      if (!selectedTrack(track1))
        continue;
      t1 = false;
      t2 = true;
      if (!selectedTrack(track2))
        continue;
      t2 = false;
      dataQA = false;

      if (produceQA) {
        registry.fill(HIST("QATrack/hSelection"), 9.5);

        registry.fill(HIST("QATrack/as/hTrack1pt"), track1.pt());
        registry.fill(HIST("QATrack/as/hTrack1DCAxy"), track1.dcaXY());
        registry.fill(HIST("QATrack/as/hTrack1DCAz"), track1.dcaZ());
        registry.fill(HIST("QATrack/as/hTrack1eta"), track1.eta());
        registry.fill(HIST("QATrack/as/hTrack1tpcNClsFound"), track1.tpcNClsFound());

        registry.fill(HIST("QATrack/as/hTrack2pt"), track2.pt());
        registry.fill(HIST("QATrack/as/hTrack2DCAxy"), track2.dcaXY());
        registry.fill(HIST("QATrack/as/hTrack2DCAz"), track2.dcaZ());
        registry.fill(HIST("QATrack/as/hTrack2eta"), track2.eta());
        registry.fill(HIST("QATrack/as/hTrack2tpcNClsFound"), track2.tpcNClsFound());
      }

      if (!selectedPair(mother, track1, track2))
        continue;

      if (produceQA) {
        registry.fill(HIST("QATrack/hSelection"), 10.5);

        registry.fill(HIST("QATrack/TPCPID/h2TracknSigma"), track1.tpcNSigmaKa(), track2.tpcNSigmaKa());
        registry.fill(HIST("QATrack/TPCPID/h2nTrack1SigmaPt"), track1.pt(), track1.tpcNSigmaKa());
        registry.fill(HIST("QATrack/TPCPID/h2nTrack2SigmaPt"), track2.pt(), track2.tpcNSigmaKa());
      }

      pointPair = fillPointPair(mother.Mag(),
                                mother.Pt(),
                                getMultiplicity(collision),
                                getCentrality(collision),
                                (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                mother.Eta(),
                                mother.Rapidity(),
                                collision.posZ(),
                                0,
                                0,
                                0);
      rsnOutput->fillUnlikepm(pointPair);
    }

    if (static_cast<bool>(produce.produceLikesign)) {

      for (const auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posDauthers, posDauthers))) {
        if (!selectedTrack(track1))
          continue;
        if (!selectedTrack(track2))
          continue;

        if (!selectedPair(mother, track1, track2))
          continue;

        if (static_cast<int>(verbose.verboselevel) > 1)
          LOGF(info, "Like-sign positive: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.Mag());

        pointPair = fillPointPair(mother.Mag(),
                                  mother.Pt(),
                                  getMultiplicity(collision),
                                  getCentrality(collision),
                                  (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                  (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                  mother.Eta(),
                                  mother.Rapidity(),
                                  collision.posZ(),
                                  0,
                                  0,
                                  0);

        rsnOutput->fillLikepp(pointPair);
      }

      for (const auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(negDauthers, negDauthers))) {
        if (!selectedTrack(track1))
          continue;
        if (!selectedTrack(track2))
          continue;

        if (!selectedPair(mother, track1, track2))
          continue;

        if (static_cast<int>(verbose.verboselevel) > 1)
          LOGF(info, "Like-sign negative: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.Mag());

        pointPair = fillPointPair(mother.Mag(),
                                  mother.Pt(),
                                  getMultiplicity(collision),
                                  getCentrality(collision),
                                  (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                  (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
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

    auto posDauthersMC = positiveMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDauthersMC = negativeMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (!collision.has_mcCollision()) {
      LOGF(warning, "No MC collision for this collision, skip...");
      return;
    }

    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersMC, negDauthersMC))) {
      registry.fill(HIST("QAMC/Truth/hMCTrack"), 0.5);

      if (!track1.has_mcParticle()) {
        LOGF(warning, "No MC particle for track, skip...");
        continue;
      }

      if (!track2.has_mcParticle()) {
        LOGF(warning, "No MC particle for track, skip...");
        continue;
      }

      MCTruthQA = true;
      t1 = true;
      if (!selectedTrack(track1))
        continue;
      t1 = false;
      t2 = true;
      if (!selectedTrack(track2))
        continue;
      t2 = false;
      MCTruthQA = false;

      registry.fill(HIST("QAMC/Truth/hMCTrack"), 9.5);

      const auto mctrack1 = track1.mcParticle();
      const auto mctrack2 = track2.mcParticle();
      int track1PDG = std::abs(mctrack1.pdgCode());
      int track2PDG = std::abs(mctrack2.pdgCode());

      if (!(track1PDG == dautherPosPDG && track2PDG == dautherNegPDG)) {
        continue;
      }
      if (produceQA)
        registry.fill(HIST("QAMC/Truth/hMCTrack"), 10.5);

      n = 0;
      for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
        for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
          if (produceQA)
            registry.fill(HIST("QAMC/Truth/hMCTrack"), 11.5);

          if (mothertrack1.pdgCode() != mothertrack2.pdgCode())
            continue;
          if (produceQA)
            registry.fill(HIST("QAMC/Truth/hMCTrack"), 12.5);

          if (mothertrack1.globalIndex() != mothertrack2.globalIndex())
            continue;
          if (produceQA)
            registry.fill(HIST("QAMC/Truth/hMCTrack"), 13.5);

          if (std::abs(mothertrack1.y()) > static_cast<float>(cut.y))
            continue;
          if (produceQA)
            registry.fill(HIST("QAMC/Truth/hMCTrack"), 14.5);

          if (std::abs(mothertrack1.pdgCode()) != motherPDG)
            continue;
          if (produceQA)
            registry.fill(HIST("QAMC/Truth/hMCTrack"), 15.5);

          if (static_cast<int>(verbose.verboselevel) > 1) {
            LOGF(info, "True: %d, d1=%d (%ld), d2=%d (%ld), mother=%d (%ld)", n, mctrack1.pdgCode(), mctrack1.globalIndex(), mctrack2.pdgCode(), mctrack2.globalIndex(), mothertrack1.pdgCode(), mothertrack1.globalIndex());
            LOGF(info, "%d px: %f, py=%f, pz=%f, px: %f, py=%f, pz=%f", n, mctrack1.px(), mctrack1.py(), mctrack1.pz(), mctrack2.px(), mctrack2.py(), mctrack2.pz());
          }

          if (!selectedPair(mother, mctrack1, mctrack2))
            continue;
          if (produceQA)
            registry.fill(HIST("QAMC/Truth/hMCTrack"), 16.5);

          if (n > 0) {
            if (produceQA)
              registry.fill(HIST("QAMC/hInvMassTrueFalse"), mother.Mag());
            continue;
          }

          pointPair = fillPointPair(mother.Mag(),
                                    mother.Pt(),
                                    getMultiplicity(collision),
                                    getCentrality(collision),
                                    (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                    (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
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
    registry.fill(HIST("QAMC/Gen/hMCEvent"), 0.5);
    if (std::abs(mcCollision.posZ()) > static_cast<float>(cut.vZ))
      return;
    registry.fill(HIST("QAMC/Gen/hMCEvent"), 1.5);

    if (collisions.size() == 0)
      return;

    for (const auto& collision : collisions) {
      registry.fill(HIST("QAMC/Gen/hMCEvent"), 2.5);

      if (!collision.has_mcCollision()) {
        LOGF(warning, "No McCollision for this collision, skip...");
        return;
      }

      auto centralityGen = 0;
      centralityGen = getCentrality(collision);
      auto multiplicityGen = 0;
      multiplicityGen = getMultiplicity(collision);

      for (const auto& particle : mcParticles) {
        registry.fill(HIST("QAMC/Gen/hMCTrack"), 0.5);

        if (std::abs(particle.y()) > static_cast<float>(cut.y))
          continue;

        registry.fill(HIST("QAMC/Gen/hMCTrack"), 1.5);

        if (particle.pdgCode() == motherPDG) {
          registry.fill(HIST("QAMC/Gen/hMCTrack"), 2.5);

          auto daughters = particle.daughters_as<aod::McParticles>();
          if (daughters.size() != 2)
            continue;

          registry.fill(HIST("QAMC/Gen/hMCTrack"), 3.5);

          auto daup = false;
          auto daun = false;

          for (const auto& dau : daughters) {
            registry.fill(HIST("QAMC/Gen/hMCTrack"), 4.5);

            if (!dau.isPhysicalPrimary())
              continue;

            registry.fill(HIST("QAMC/Gen/hMCTrack"), 5.5);

            if (dau.pdgCode() == dautherPosPDG) {
              daup = true;
              d1.SetXYZM(dau.px(), dau.py(), dau.pz(), massPos);
            } else if (dau.pdgCode() == -dautherNegPDG) {
              daun = true;
              d2.SetXYZM(dau.px(), dau.py(), dau.pz(), massNeg);
            }
          }
          if (!daup && !daun)
            continue;

          registry.fill(HIST("QAMC/Gen/hMCTrack"), 6.5);

          mother = d1 + d2;

          pointPair = fillPointPair(mother.Mag(),
                                    mother.Pt(),
                                    multiplicityGen,
                                    centralityGen,
                                    std::abs(static_cast<float>(cut.tpcnSigmaPos) / 2.0),
                                    std::abs(static_cast<float>(cut.tpcnSigmaNeg) / 2.0),
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
  PROCESS_SWITCH(PhianalysisTHnSparse, processGen, "Process MC Mateched.", true);

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

        auto posDauthersc1 = positive->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto posDauthersc2 = positive->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);
        auto negDauthersc1 = negative->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto negDauthersc2 = negative->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);

        if (produceQA) {
          registry.fill(HIST("QAMixing/h2mu1_mu2"), getMultiplicity(c1), getMultiplicity(c2));
          registry.fill(HIST("QAMixing/h2ce1_ce2"), getCentrality(c1), getCentrality(c2));
          registry.fill(HIST("QAMixing/h2vz1_vz2"), c1.posZ(), c2.posZ());
        }

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersc1, negDauthersc2))) {
          if (produceQA)
            registry.fill(HIST("QAMixing/hTrackSelection"), 0.5);

          if (!selectedTrack(track1))
            continue;
          if (produceQA)
            registry.fill(HIST("QAMixing/hTrackSelection"), 1.5);
          if (!selectedTrack(track2))
            continue;
          if (produceQA)
            registry.fill(HIST("QAMixing/hTrackSelection"), 2.5);

          if (!selectedPair(mother, track1, track2))
            continue;
          if (produceQA)
            registry.fill(HIST("QAMixing/hTrackSelection"), 3.5);

          pointPair = fillPointPair(mother.Mag(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                    (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    c1.posZ(),
                                    getMultiplicity(c2),
                                    getCentrality(c2),
                                    c2.posZ());

          rsnOutput->fillMixingpm(pointPair);
        }

        if (static_cast<bool>(produce.produceLikesign)) {

          for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersc1, posDauthersc2))) {

            if (!selectedTrack(track1))

              continue;
            if (!selectedTrack(track2))
              continue;

            if (!selectedPair(mother, track1, track2))
              continue;

            pointPair = fillPointPair(mother.Mag(),
                                      mother.Pt(),
                                      getMultiplicity(c1),
                                      getCentrality(c1),
                                      (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                      (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                      mother.Eta(),
                                      mother.Rapidity(),
                                      c1.posZ(),
                                      getMultiplicity(c2),
                                      getCentrality(c2),
                                      c2.posZ());

            rsnOutput->fillMixingpp(pointPair);
          }

          for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(negDauthersc1, negDauthersc2))) {

            if (!selectedTrack(track1))

              continue;
            if (!selectedTrack(track2))
              continue;

            if (!selectedPair(mother, track1, track2))
              continue;
            pointPair = fillPointPair(mother.Mag(),
                                      mother.Pt(),
                                      getMultiplicity(c1),
                                      getCentrality(c1),
                                      (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                      (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                      mother.Eta(),
                                      mother.Rapidity(),
                                      c1.posZ(),
                                      getMultiplicity(c2),
                                      getCentrality(c2),
                                      c2.posZ());

            rsnOutput->fillMixingmm(pointPair);
          }
        }

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersc2, negDauthersc1))) {

          if (!selectedTrack(track1))

            continue;
          if (!selectedTrack(track2))
            continue;

          if (!selectedPair(mother, track1, track2))
            continue;

          pointPair = fillPointPair(mother.Mag(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                    (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
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

        auto posDauthersc1 = positive->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto posDauthersc2 = positive->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);
        auto negDauthersc1 = negative->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
        auto negDauthersc2 = negative->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);

        if (produceQA) {
          registry.fill(HIST("QAMixing/h2mu1_mu2"), getMultiplicity(c1), getMultiplicity(c2));
          registry.fill(HIST("QAMixing/h2ce1_ce2"), getCentrality(c1), getCentrality(c2));
          registry.fill(HIST("QAMixing/h2vz1_vz2"), c1.posZ(), c2.posZ());
        }

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersc1, negDauthersc2))) {
          if (produceQA)
            registry.fill(HIST("QAMixing/hTrackSelection"), 0.5);

          if (!selectedTrack(track1))
            continue;
          if (produceQA)
            registry.fill(HIST("QAMixing/hTrackSelection"), 1.5);
          if (!selectedTrack(track2))
            continue;
          if (produceQA)
            registry.fill(HIST("QAMixing/hTrackSelection"), 2.5);

          if (!selectedPair(mother, track1, track2))
            continue;
          if (produceQA)
            registry.fill(HIST("QAMixing/hTrackSelection"), 3.5);

          pointPair = fillPointPair(mother.Mag(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                    (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    c1.posZ(),
                                    getMultiplicity(c2),
                                    getCentrality(c2),
                                    c2.posZ());

          rsnOutput->fillMixingpm(pointPair);
        }

        if (static_cast<bool>(produce.produceLikesign)) {

          for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersc1, posDauthersc2))) {

            if (!selectedTrack(track1))

              continue;
            if (!selectedTrack(track2))
              continue;

            if (!selectedPair(mother, track1, track2))
              continue;

            pointPair = fillPointPair(mother.Mag(),
                                      mother.Pt(),
                                      getMultiplicity(c1),
                                      getCentrality(c1),
                                      (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                      (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                      mother.Eta(),
                                      mother.Rapidity(),
                                      c1.posZ(),
                                      getMultiplicity(c2),
                                      getCentrality(c2),
                                      c2.posZ());

            rsnOutput->fillMixingpp(pointPair);
          }

          for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(negDauthersc1, negDauthersc2))) {

            if (!selectedTrack(track1))

              continue;
            if (!selectedTrack(track2))
              continue;

            if (!selectedPair(mother, track1, track2))
              continue;
            pointPair = fillPointPair(mother.Mag(),
                                      mother.Pt(),
                                      getMultiplicity(c1),
                                      getCentrality(c1),
                                      (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                      (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                      mother.Eta(),
                                      mother.Rapidity(),
                                      c1.posZ(),
                                      getMultiplicity(c2),
                                      getCentrality(c2),
                                      c2.posZ());

            rsnOutput->fillMixingmm(pointPair);
          }
        }

        for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersc2, negDauthersc1))) {

          if (!selectedTrack(track1))

            continue;
          if (!selectedTrack(track2))
            continue;

          if (!selectedPair(mother, track1, track2))
            continue;

          pointPair = fillPointPair(mother.Mag(),
                                    mother.Pt(),
                                    getMultiplicity(c1),
                                    getCentrality(c1),
                                    (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                    (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
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
  PROCESS_SWITCH(PhianalysisTHnSparse, processMixed, "Process Mixing Event.", true);

  void processGenOld(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    if (!static_cast<bool>(produce.produceTrue))
      return;

    registry.fill(HIST("QAMC/hMC"), 0.5);

    if (std::abs(mcCollision.posZ()) > static_cast<float>(cut.vZ))
      return;

    registry.fill(HIST("QAMC/hMC"), 1.5);

    for (const auto& particle : mcParticles) {
      registry.fill(HIST("QAMC/hMC"), 2.5);
      if (std::abs(particle.y()) > static_cast<float>(cut.y))
        continue;

      registry.fill(HIST("QAMC/hMC"), 3.5);

      if (particle.pdgCode() == motherPDG) {
        auto daughters = particle.daughters_as<aod::McParticles>();
        if (daughters.size() != 2)
          continue;

        registry.fill(HIST("QAMC/hMC"), 4.5);

        auto daup = false;
        auto daun = false;

        for (const auto& dau : daughters) {
          if (!dau.isPhysicalPrimary())
            continue;

          if (dau.pdgCode() == dautherPosPDG) {
            daup = true;
            d1.SetXYZM(dau.px(), dau.py(), dau.pz(), massPos);
          } else if (dau.pdgCode() == -dautherNegPDG) {
            daun = true;
            d2.SetXYZM(dau.px(), dau.py(), dau.pz(), massNeg);
          }
        }
        if (!daup && !daun)
          continue;

        registry.fill(HIST("QAMC/hMC"), 5.5);

        mother = d1 + d2;

        pointPair = fillPointPair(mother.Mag(),
                                  mother.Pt(),
                                  0,
                                  0,
                                  std::abs(static_cast<float>(cut.tpcnSigmaPos) / 2.0),
                                  std::abs(static_cast<float>(cut.tpcnSigmaNeg) / 2.0),
                                  mother.Eta(),
                                  mother.Rapidity(),
                                  mcCollision.posZ(),
                                  0,
                                  0,
                                  0);

        rsnOutput->fillUnlikegenOld(pointPair);
      }
    }
  }
  PROCESS_SWITCH(PhianalysisTHnSparse, processGenOld, "Process generated.", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PhianalysisTHnSparse>(cfgc)};
}
