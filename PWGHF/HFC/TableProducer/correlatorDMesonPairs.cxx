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

/// \file correlatorDMesonPairs.cxx
/// \brief D0(bar) correlator task - data-like, MC-reco and MC-kine analyses.
///
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseD0ToKPi.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFC/DataModel/DMesonPairsTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <THnSparse.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
using DecayChannel = o2::hf_decay::hf_cand_2prong::DecayChannelMain;

enum CandidateType {
  SelectedD = 0, // This particle is selected as a D
  SelectedDbar,  // This particle is selected as a Dbar
  TrueD,         // This particle is a true D
  TrueDbar       // This particle is a true Dbar
};

enum PairTypeOfSelMassSel {
  DD = 0,   // This is a D0-D0 pair
  DbarDbar, // This is a D0bar-D0bar pair
  DDbar,    // This is a D0-D0bar pair
  DbarD     // This is a D0bar-D0 pair
};
} // namespace

using McParticlesPlus2Prong = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;
// definition of ME variables
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;

struct HfCorrelatorD0HadronsSelection {
  Produces<aod::DmesonSelection> d0Sel;

  Configurable<bool> useSel8{"useSel8", true, "Flag for applying sel8 for collision selection"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<bool> doSelD0Collision{"doSelD0Collision", true, "Select collisions with at least one D0"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  SliceCache cache;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CandidatesD0Data = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;

  Partition<CandidatesD0Data> selectedD0Candidates =
    aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 ||
    aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  // Process function to select collisions with at least one D0 candidate passing the selection criteria
  void processD0SelectionData(SelCollisions::iterator const& collision,
                              CandidatesD0Data const&)
  {
    bool const isSel8 = !useSel8 || collision.sel8();
    bool const isNoSameBunchPileUp = !selNoSameBunchPileUpColl ||
                                     collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
    bool isD0Found = !doSelD0Collision;

    if (doSelD0Collision) {
      auto grouped = selectedD0Candidates.sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      for (const auto& candidate : grouped) {
        if (std::abs(HfHelper::yD0(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          continue;
        }
        isD0Found = true;
        break;
      }
    }

    d0Sel(isD0Found && isSel8 && isNoSameBunchPileUp);
  }
  PROCESS_SWITCH(HfCorrelatorD0HadronsSelection, processD0SelectionData, "Process D0 Selection Data", true);
};

struct HfCorrelatorDMesonPairs {

  Produces<aod::D0Pair> entryD0Pair;
  Produces<aod::D0PairMl> entryD0PairMl;
  Produces<aod::D0PairMcInfo> entryD0PairMcInfo;
  Produces<aod::D0PairMcGen> entryD0PairMcGen;
  Produces<aod::D0PairMcGenInfo> entryD0PairMcGenInfo;

  // Tables for event mixing
  Produces<aod::DMesonInfoCand1> entryDMesonInfoCand1;
  Produces<aod::DMesonInfoCand2> entryDMesonInfoCand2;
  Produces<aod::AssocHadInfo> entryAssocHad;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "maximum |y| of D0 candidates"};
  Configurable<float> ptCandMin{"ptCandMin", -1., "minimum pT of D0 candidates"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots"};
  Configurable<bool> selectSignalRegionOnly{"selectSignalRegionOnly", false, "only use events close to PDG peak"};
  Configurable<float> massCut{"massCut", 0.05, "Maximum deviation from PDG peak allowed for signal region"};
  Configurable<bool> removeAmbiguous{"removeAmbiguous", false, "Flag to remove ambiguous candidates"};
  Configurable<float> ptMaxRemoveAmbiguous{"ptMaxRemoveAmbiguous", 5.0, "Max. pT to remove the ambiguous candidates"};

  // Event mixing configurables
  Configurable<bool> applyMixedEvent{"applyMixedEvent", false, "Flag to make the event mixing"};
  Configurable<bool> daughterTracksCutFlag{"daughterTracksCutFlag", false, "Flag to add cut on daughter tracks"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCA_xy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCA_z of tracks"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<int> multMin{"multMin", 0, "minimum multiplicity accepted"};
  Configurable<int> multMax{"multMax", 10000, "maximum multiplicity accepted"};

  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {static_cast<const double*>(hf_cuts_ml::Cuts[0]), hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};

  // ML model CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTD0"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_D0ToKPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  // Table aliases. The partitions below must have exactly the same type as the
  // candidate table taken by the process functions, otherwise they are never bound.
  using SelCollisionsWithD0 = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::DmesonSelection>>;
  using CandidatesD0Data = soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0>;
  using CandidatesD0McRec = soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfCand2ProngMcRec>;
  using TracksData = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection>>;

  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);

  SliceCache cache;

  Partition<CandidatesD0Data> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<CandidatesD0McRec> selectedD0CandidatesMc = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  o2::analysis::HfMlResponseD0ToKPi<float> hfMlResponse;
  o2::ccdb::CcdbApi ccdbApi;

  std::vector<float> outputMlD0Cand1;
  std::vector<float> outputMlD0barCand1;

  std::vector<float> outputMlD0Cand2;
  std::vector<float> outputMlD0barCand2;

  // ThnSparse for ML outputScores and Vars
  ConfigurableAxis thnConfigAxisBkgScore{"thnConfigAxisBkgScore", {100, 0, 1}, "Bkg score bins"};
  ConfigurableAxis thnConfigAxisSignalScore{"thnConfigAxisSignalScore", {100, 0, 1}, "Signal score bins"};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {120, 1.5848, 2.1848}, "Cand. inv-mass bins"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {500, 0, 50}, "Cand. pT bins"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, "Cand. origin type"};
  ConfigurableAxis thnConfigAxisCandType{"thnConfigAxisCandType", {6, -0.5, 5.5}, "D0 type"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};

  HistogramConfigSpec hTH1Pt{HistType::kTH1F, {{180, 0., 36.}}};
  HistogramConfigSpec hTH1Y{HistType::kTH1F, {{100, -5., 5.}}};
  HistogramConfigSpec hTH1NContrib{HistType::kTH1F, {{200, -0.5, 199.5}}};
  HistogramConfigSpec hTH1Phi{HistType::kTH1F, {{32, 0., o2::constants::math::TwoPI}}};
  HistogramConfigSpec hTH2Pid{HistType::kTH2F, {{500, 0., 10.}, {400, -20., 20.}}};
  HistogramConfigSpec hTH3PtVsYVsNContrib{HistType::kTH3F, {{360, 0., 36.}, {20, -1., 1.}, {120, -0.5, 119.5}}};

  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 100000.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis binsZVtx{"binsZVtx", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "Mixing bins - z-vertex"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};
  std::unique_ptr<BinningType> corrBinning;

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "D meson candidates;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtCandAfterCut", "D meson candidates after pT cut;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtProng0", "D meson candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtProng1", "D meson candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hEta", "D meson candidates;candidate #it{#eta};entries", hTH1Y},
     {"hPhi", "D meson candidates;candidate #it{#varphi};entries", hTH1Phi},
     {"hY", "D meson candidates;candidate #it{y};entries", hTH1Y},
     {"hPVContrib", "D meson candidates;candidate Number of PV contributors;entries", hTH1NContrib},
     // MC Gen plots
     {"hPtCandMcGen", "D meson candidates MC Gen;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hPtCandAfterCutMcGen", "D meson candidates after pT cut;candidate #it{p}_{T} (GeV/#it{c});entries", hTH1Pt},
     {"hEtaMcGen", "D meson candidates MC Gen;candidate #it{#eta};entries", hTH1Y},
     {"hPhiMcGen", "D meson candidates MC Gen;candidate #it{#varphi};entries", hTH1Phi},
     {"hPtVsYVsNContribMcGen", "D meson candidates MC Gen;candidate #it{p}_{T} (GeV/#it{c});#it{y};Number of contributors", hTH3PtVsYVsNContrib},
     {"hPtVsYVsNContribMcGenPrompt", "D meson candidates MC Gen Prompt;candidate #it{p}_{T} (GeV/#it{c});#it{y};Number of contributors", hTH3PtVsYVsNContrib},
     {"hPtVsYVsNContribMcGenNonPrompt", "D meson candidates MC Gen Non-prompt;candidate #it{p}_{T} (GeV/#it{c});#it{y};Number of contributors", hTH3PtVsYVsNContrib},
     {"hNContribMcGen", "D meson candidates MC Gen;Number of PV contributors", hTH1NContrib},
     // MC Rec plots
     {"hPtVsYVsNContribMcRec", "D meson candidates MC Rec;candidate #it{p}_{T} (GeV/#it{c});#it{y};Number of contributors", hTH3PtVsYVsNContrib},
     {"hPtVsYVsNContribMcRecPrompt", "D meson candidates MC Rec Prompt;candidate #it{p}_{T} (GeV/#it{c});#it{y};Number of contributors", hTH3PtVsYVsNContrib},
     {"hPtVsYVsNContribMcRecNonPrompt", "D meson candidates MC Rec Non-prompt;candidate #it{p}_{T} (GeV/#it{c});#it{y};Number of contributors", hTH3PtVsYVsNContrib},
     {"hNContribMcRec", "D meson candidates MC Rec;Number of PV contributors", hTH1NContrib},
     // PID plots ----- Not definitively here
     {"PID/hTofNSigmaPi", "(TOFsignal-time#pi)/tofSigPid;p[GeV/c];(TOFsignal-time#pi)/tofSigPid", hTH2Pid},
     {"PID/hTofNSigmaKa", "(TOFsignal-timeK)/tofSigPid;p[GeV/c];(TOFsignal-timeK)/tofSigPid", hTH2Pid},
     {"PID/hTpcNSigmaPi", "(TPCsignal-time#pi)/tpcSigPid;p[GeV/c];(TPCsignal-time#pi)/tpcSigPid", hTH2Pid},
     {"PID/hTpcNSigmaKa", "(TPCsignal-timeK)/tpcSigPid;p[GeV/c];(TPCsignal-timeK)/tpcSigPid", hTH2Pid},
     {"PID/hTpcTofNSigmaPi", "(TPC+TOFsignal-time#pi)/tpcTofSigPid;p[GeV/#it{c}];(TPC+TOFsignal-time#pi)/tpcTofSigPid", hTH2Pid},
     {"PID/hTpcTofNSigmaKa", "(TPC+TOFsignal-timeK)/tpcTofSigPid;p[GeV/c];(TPC+TOFsignal-timeK)/tpcTofSigPid", hTH2Pid}}};

  void init(InitContext&)
  {
    if (applyMl) {
      hfMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      hfMlResponse.init();
    }

    auto vbins = static_cast<std::vector<double>>(binsPt);
    constexpr int NBinsSelStatus = 25;
    std::array<std::string, NBinsSelStatus> labels;

    labels[0] = "total # of Selected pairs";
    // Cand1 analysis
    labels[1] = "total # of Selected Cands 1";
    labels[2] = "# of Selected D Cand 1 ONLY";
    labels[3] = "# of Selected Dbar Cand 1 ONLY";
    labels[4] = "# of Selected Simultaneous D + Dbar Cand 1";
    labels[5] = "# of True D Cand 1";
    labels[6] = "# of True Dbar Cand 1";
    // Cand2 analysis
    labels[7] = "total # of Selected Cands 2";
    labels[8] = "# of Selected D Cand 2 ONLY";
    labels[9] = "# of Selected Dbar Cand 2 ONLY";
    labels[10] = "# of Selected Simultaneous D + Dbar Cand 2";
    labels[11] = "# of True D Cand 2";
    labels[12] = "# of True Dbar Cand 2";
    // Pair analysis
    labels[13] = "# of D+D Pairs";
    labels[14] = "# of Dbar+Dbar Pairs";
    labels[15] = "# of D+Dbar Pairs";
    labels[16] = "# of Dbar+D Pairs";
    labels[17] = "# of D+D ONLY Pairs";
    labels[18] = "# of Dbar+Dbar ONLY Pairs";
    labels[19] = "# of D+Dbar ONLY Pairs";
    labels[20] = "# of Dbar+D ONLY Pairs";
    // True pair analysis
    labels[21] = "# of True D+D Pairs";
    labels[22] = "# of True Dbar+Dbar Pairs";
    labels[23] = "# of True D+Dbar Pairs";
    labels[24] = "# of True Dbar+D Pairs";

    AxisSpec const axisSelStatus = {NBinsSelStatus, 0.5, NBinsSelStatus + 0.5, ""};
    registry.add("hSelectionStatus", "D Meson candidates;selection status;entries", HistType::kTH1F, {axisSelStatus});
    registry.add("hSelectionStatusMcGen", "D Meson candidates MC Gen;selection status;entries", HistType::kTH1F, {axisSelStatus});

    for (int iBin = 0; iBin < NBinsSelStatus; iBin++) {
      registry.get<TH1>(HIST("hSelectionStatus"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      registry.get<TH1>(HIST("hSelectionStatusMcGen"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    constexpr int NBinsMatching = 8;
    std::array<std::string, NBinsMatching> labelsMatching;
    // Cand1 analysis
    labelsMatching[0] = "total # of Cand 1";
    labelsMatching[1] = "# of matched D Cand 1";
    labelsMatching[2] = "# of matched Dbar Cand 1";
    labelsMatching[3] = "# of unmatched Cand 1";
    // Cand2 analysis
    labelsMatching[4] = "total # of Cand 2";
    labelsMatching[5] = "# of matched D Cand 2";
    labelsMatching[6] = "# of matched Dbar Cand 2";
    labelsMatching[7] = "# of unmatched Cand 2";

    AxisSpec const axisMatching = {NBinsMatching, 0.5, NBinsMatching + 0.5, ""};
    registry.add("hMatchingMcRec", "D Meson candidates; MC matching status;entries", HistType::kTH1F, {axisMatching});
    registry.add("hMatchingMcGen", "D Meson candidates; MC matching status;entries", HistType::kTH1F, {axisMatching});

    for (int iBin = 0; iBin < NBinsMatching; iBin++) {
      registry.get<TH1>(HIST("hMatchingMcRec"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMatching[iBin].data());
      registry.get<TH1>(HIST("hMatchingMcGen"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMatching[iBin].data());
    }

    constexpr int NBinsSinglePart = 6;
    std::array<std::string, NBinsSinglePart> labelsSinglePart;
    // Candidate analysis
    labelsSinglePart[0] = "total # of Candidates";
    labelsSinglePart[1] = "# of selected D";
    labelsSinglePart[2] = "# of selected Dbar";
    labelsSinglePart[3] = "# of selected D and Dbar";
    labelsSinglePart[4] = "# of true D";
    labelsSinglePart[5] = "# of true Dbar";

    AxisSpec const axisSinglePart = {NBinsSinglePart, 0.5, NBinsSinglePart + 0.5, ""};
    registry.add("hStatusSinglePart", "D Meson candidates; MC matching status;entries", HistType::kTH1F, {axisSinglePart});
    registry.add("hStatusSinglePartMcGen", "D Meson candidates; MC matching status;entries", HistType::kTH1F, {axisSinglePart});

    for (int iBin = 0; iBin < NBinsSinglePart; iBin++) {
      registry.get<TH1>(HIST("hStatusSinglePart"))->GetXaxis()->SetBinLabel(iBin + 1, labelsSinglePart[iBin].data());
      registry.get<TH1>(HIST("hStatusSinglePartMcGen"))->GetXaxis()->SetBinLabel(iBin + 1, labelsSinglePart[iBin].data());
    }

    AxisSpec const axisInputD0 = {200, -0.5, 199.5};
    registry.add("hInputCheckD0", "Check on input D0 meson candidates/event", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0bar", "Check on input D0bar meson candidates/event", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0AndD0bar", "Check on input D0 & D0bar meson candidates/event", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0OrD0bar", "Check on input D0 | D0bar meson candidates/event", {HistType::kTH1F, {axisInputD0}});
    // MC Gen
    registry.add("hInputCheckD0McGen", "Check on input D0 meson candidates/event MC Gen", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0barMcGen", "Check on input D0bar meson candidates/event MC Gen", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0AndD0barMcGen", "Check on input D0 & D0bar meson candidates/event MC Gen", {HistType::kTH1F, {axisInputD0}});
    registry.add("hInputCheckD0OrD0barMcGen", "Check on input D0 | D0bar meson candidates/event MC Gen", {HistType::kTH1F, {axisInputD0}});

    registry.add("hMass", "D Meson pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassSig", "MC-matched D Meson pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassMcRecPrompt", "D Meson pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassMcRecNonPrompt", "D Meson pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassMcRecReflections", "D Meson pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (#pi K) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisY{thnConfigAxisY, "y"};
    const AxisSpec thnAxisOrigin{thnConfigAxisOrigin, "Origin"};
    const AxisSpec thnAxisCandType{thnConfigAxisCandType, "D0 type"};
    const AxisSpec thnAxisNumPvContr{thnConfigAxisNumPvContr, "Number of PV contributors"};

    std::vector<AxisSpec> axes = {thnAxisMass, thnAxisPt, thnAxisY, thnAxisNumPvContr, thnAxisOrigin, thnAxisCandType};
    if (applyMl) {
      const AxisSpec thnAxisBkgScore{thnConfigAxisBkgScore, "BDT score bkg"};
      const AxisSpec thnAxisSignalScore{thnConfigAxisSignalScore, "BDT score signal"};

      axes.insert(axes.begin(), thnAxisSignalScore);
      axes.insert(axes.begin(), thnAxisBkgScore);

      registry.add("hnDMesonMl", "THn for D Meson candidates", HistType::kTHnSparseD, axes);
      registry.get<THnSparse>(HIST("hnDMesonMl"))->Sumw2();
    } else {
      registry.add("hnDMeson", "Thn for D0 candidates", HistType::kTHnSparseD, axes);
      registry.get<THnSparse>(HIST("hnDMeson"))->Sumw2();
    }

    // Event mixing
    AxisSpec const axisPoolBin = {binsPoolBin, "PoolBin"};
    corrBinning = std::make_unique<BinningType>(BinningType{{binsZVtx, binsMultiplicity}});

    if (applyMixedEvent) {
      registry.add("hMultFT0M", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
      registry.add("hZvtx", "z vertex;z vertex;entries", {HistType::kTH1F, {{200, -20., 20.}}});
      registry.add("hD0Bin", "D0 selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{axisPoolBin}}});
      registry.add("hD0Cand2Bin", "D0 selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{axisPoolBin}}});
    }
  }

  /// Kinematic, signal-region and ambiguity selections of a candidate
  /// \return true if the candidate passes the selections
  template <typename T>
  bool passKinematicSelections(const T& candidate, const bool isDCand, const bool isDbarCand)
  {
    if (std::abs(HfHelper::yD0(candidate)) > yCandMax) {
      return false;
    }
    if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
      return false;
    }
    bool const isSignalD0Cand = std::abs(HfHelper::invMassD0ToPiK(candidate) - MassD0) < massCut;
    bool const isSignalD0barCand = std::abs(HfHelper::invMassD0barToKPi(candidate) - MassD0Bar) < massCut;
    if (selectSignalRegionOnly && !(isSignalD0Cand || isSignalD0barCand)) {
      return false;
    }
    // Remove ambiguous D0 candidates if flag is true
    if (removeAmbiguous && (isDCand && isDbarCand) && candidate.pt() < ptMaxRemoveAmbiguous) {
      return false;
    }
    return true;
  }

  /// Sets bits to select candidate type for D0
  /// SelectedD and SelectedDbar bits look at whether the candidate passed the selection flags.
  /// \param candidate is candidate
  /// \return bitmap with type of candidate
  template <bool IsMcRec, typename T>
  uint8_t assignCandidateTypeD0(const T& candidate)
  {
    uint8_t candidateType(0);
    if (candidate.isSelD0() >= selectionFlagD0) {
      SETBIT(candidateType, SelectedD);
    }
    if (candidate.isSelD0bar() >= selectionFlagD0bar) {
      SETBIT(candidateType, SelectedDbar);
    }
    if constexpr (IsMcRec) {
      if (candidate.flagMcMatchRec() == DecayChannel::D0ToPiK) { // matched as D0
        SETBIT(candidateType, TrueD);
      }
      if (candidate.flagMcMatchRec() == -DecayChannel::D0ToPiK) { // matched as D0bar
        SETBIT(candidateType, TrueDbar);
      }
    }
    return candidateType;
  }

  /// Sets bits to select candidate type at generator level
  /// \param particle is particle
  /// \return bitmap with type of gen particle - they are selected as True
  template <typename T>
  uint8_t assignParticleTypeD0Gen(const T& particle)
  {
    uint8_t particleType(0);
    if (particle.pdgCode() == Pdg::kD0) { // just checking the particle PDG, not the decay channel
      SETBIT(particleType, TrueD);
    }
    if (particle.pdgCode() == (-Pdg::kD0)) { // just checking the particle PDG, not the decay channel
      SETBIT(particleType, TrueDbar);
    }
    return particleType;
  }

  /// Check if the candidate is a D
  bool isD(const uint8_t candidateType)
  {
    return (TESTBIT(candidateType, SelectedD));
  }

  /// Check if the candidate is a Dbar
  bool isDbar(const uint8_t candidateType)
  {
    return (TESTBIT(candidateType, SelectedDbar));
  }

  /// Check if the candidate is a true D
  bool isTrueD(const uint8_t candidateType)
  {
    return (TESTBIT(candidateType, TrueD));
  }

  /// Check if the candidate is a true Dbar
  bool isTrueDbar(const uint8_t candidateType)
  {
    return (TESTBIT(candidateType, TrueDbar));
  }

  /// Rapidity and pT selections of gen-level particles
  template <typename T>
  bool passGenSelections(const T& particle)
  {
    return std::abs(particle.y()) <= yCandMax && (ptCandMin < 0. || particle.pt() >= ptCandMin);
  }

  /// Canonical pair ordering: D0 before D0bar; for equal flavours, higher pT first.
  /// \return true if candidate 1 and candidate 2 should be exchanged
  template <typename T>
  static bool shouldSwapOrder(const bool isD1, const bool isDbar1, const bool isD2, const bool isDbar2, const T& candidate1, const T& candidate2)
  {
    if (isD1 != isD2 || isDbar1 != isDbar2) { // different flavour content
      return !isD1 && isD2;                    // Dbar(-only) first, D second -> swap
    }
    return candidate1.pt() < candidate2.pt();
  }

  /// ML selections of a candidate
  /// \return true if the candidate is selected by the ML model (as D0 and/or D0bar)
  template <typename T>
  bool passMlSelections(const T& candidate, const bool isDCand, const bool isDbarCand, std::vector<float>& outputMlD0Cand, std::vector<float>& outputMlD0barCand)
  {
    bool isSelectedMlD0Cand = false;
    bool isSelectedMlD0barCand = false;

    if (isDCand) {
      std::vector<float> inputFeaturesD0 = hfMlResponse.getInputFeatures<false>(candidate, o2::constants::physics::kD0);
      isSelectedMlD0Cand = hfMlResponse.isSelectedMl(inputFeaturesD0, candidate.pt(), outputMlD0Cand);
    }
    if (isDbarCand) {
      std::vector<float> inputFeaturesD0bar = hfMlResponse.getInputFeatures<false>(candidate, o2::constants::physics::kD0Bar);
      isSelectedMlD0barCand = hfMlResponse.isSelectedMl(inputFeaturesD0bar, candidate.pt(), outputMlD0barCand);
    }

    return isSelectedMlD0Cand || isSelectedMlD0barCand;
  }

  /// Fill PID related plots for D0 and D0bar
  /// \param candidate is candidate
  template <typename T>
  void analysePid(const T& candidate)
  {
    if (candidate.isSelD0() >= selectionFlagD0) {
      registry.fill(HIST("PID/hTofNSigmaPi"), candidate.ptProng0(), candidate.nSigTofPi0());
      registry.fill(HIST("PID/hTofNSigmaKa"), candidate.ptProng1(), candidate.nSigTofKa1());
      registry.fill(HIST("PID/hTpcNSigmaPi"), candidate.ptProng0(), candidate.nSigTpcPi0());
      registry.fill(HIST("PID/hTpcNSigmaKa"), candidate.ptProng1(), candidate.nSigTpcKa1());
      registry.fill(HIST("PID/hTpcTofNSigmaPi"), candidate.ptProng0(), candidate.tpcTofNSigmaPi0());
      registry.fill(HIST("PID/hTpcTofNSigmaKa"), candidate.ptProng1(), candidate.tpcTofNSigmaKa1());
    }
    if (candidate.isSelD0bar() >= selectionFlagD0bar) {
      registry.fill(HIST("PID/hTofNSigmaPi"), candidate.ptProng1(), candidate.nSigTofPi1());
      registry.fill(HIST("PID/hTofNSigmaKa"), candidate.ptProng0(), candidate.nSigTofKa0());
      registry.fill(HIST("PID/hTpcNSigmaPi"), candidate.ptProng1(), candidate.nSigTpcPi1());
      registry.fill(HIST("PID/hTpcNSigmaKa"), candidate.ptProng0(), candidate.nSigTpcKa0());
      registry.fill(HIST("PID/hTpcTofNSigmaPi"), candidate.ptProng1(), candidate.tpcTofNSigmaPi1());
      registry.fill(HIST("PID/hTpcTofNSigmaKa"), candidate.ptProng0(), candidate.tpcTofNSigmaKa0());
    }
  }

  /// Fill counters for D0 and D0bar
  /// \param d0Candidates contains the D0 candidates of the collision
  template <bool IsMcRec, typename T>
  void getCountersPerEvent(const T& d0Candidates)
  {
    int nDevent = 0, nDbarevent = 0, nDDbarevent = 0, nDorDbarevent = 0;
    for (const auto& candidate : d0Candidates) {
      // Get counters per event
      auto const candidateType = assignCandidateTypeD0<IsMcRec>(candidate); // Candidate type attribution
      bool const isDCand = isD(candidateType);
      bool const isDbarCand = isDbar(candidateType);
      bool const isTrueDCand = isTrueD(candidateType);
      bool const isTrueDbarCand = isTrueDbar(candidateType);

      if (!passKinematicSelections(candidate, isDCand, isDbarCand)) {
        continue;
      }

      if (isDCand) {
        nDevent++;
      }
      if (isDbarCand) {
        nDbarevent++;
      }
      if (isDCand && isDbarCand) {
        nDDbarevent++;
      }
      if (isDCand || isDbarCand) {
        nDorDbarevent++;
      }

      // Fill selection status single particle
      registry.fill(HIST("hStatusSinglePart"), 1);
      if (isDCand && !isDbarCand) {
        registry.fill(HIST("hStatusSinglePart"), 2);
      } else if (isDbarCand && !isDCand) {
        registry.fill(HIST("hStatusSinglePart"), 3);
      } else if (isDCand && isDbarCand) {
        registry.fill(HIST("hStatusSinglePart"), 4);
      }
      // fill true info (available only in MC)
      if (isTrueDCand) {
        registry.fill(HIST("hStatusSinglePart"), 5);
      } else if (isTrueDbarCand) {
        registry.fill(HIST("hStatusSinglePart"), 6);
      }
    }
    if (nDevent > 0) {
      registry.fill(HIST("hInputCheckD0"), nDevent);
    }
    if (nDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0bar"), nDbarevent);
    }
    if (nDDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0AndD0bar"), nDDbarevent);
    }
    if (nDorDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0OrD0bar"), nDorDbarevent);
    }
  }

  void fillSelectionStatusHisto(const uint8_t candidateType1, const uint8_t candidateType2)
  {
    bool const isDCand1 = isD(candidateType1), isDCand2 = isD(candidateType2);
    bool const isDbarCand1 = isDbar(candidateType1), isDbarCand2 = isDbar(candidateType2);
    bool const isTrueDCand1 = isTrueD(candidateType1), isTrueDCand2 = isTrueD(candidateType2);
    bool const isTrueDbarCand1 = isTrueDbar(candidateType1), isTrueDbarCand2 = isTrueDbar(candidateType2);

    /// Fill information on the D candidates
    registry.fill(HIST("hSelectionStatus"), 2);
    if (isDCand1 && !isDbarCand1) {
      registry.fill(HIST("hSelectionStatus"), 3);
    } else if (isDbarCand1 && !isDCand1) {
      registry.fill(HIST("hSelectionStatus"), 4);
    } else if (isDCand1 && isDbarCand1) {
      registry.fill(HIST("hSelectionStatus"), 5);
    }

    registry.fill(HIST("hSelectionStatus"), 8);
    if (isDCand2 && !isDbarCand2) {
      registry.fill(HIST("hSelectionStatus"), 9);
    } else if (isDbarCand2 && !isDCand2) {
      registry.fill(HIST("hSelectionStatus"), 10);
    } else if (isDCand2 && isDbarCand2) {
      registry.fill(HIST("hSelectionStatus"), 11);
    }

    /// Collect information on the D pairs
    registry.fill(HIST("hSelectionStatus"), 1);
    if (isDCand1 && isDCand2) {
      registry.fill(HIST("hSelectionStatus"), 14);
      if ((!isDbarCand1) && (!isDbarCand2)) {
        registry.fill(HIST("hSelectionStatus"), 18);
      }
    }
    if (isDbarCand1 && isDbarCand2) {
      registry.fill(HIST("hSelectionStatus"), 15);
      if ((!isDCand1) && (!isDCand2)) {
        registry.fill(HIST("hSelectionStatus"), 19);
      }
    }
    if (isDCand1 && isDbarCand2) {
      registry.fill(HIST("hSelectionStatus"), 16);
      if ((!isDbarCand1) && (!isDCand2)) {
        registry.fill(HIST("hSelectionStatus"), 20);
      }
    }
    if (isDbarCand1 && isDCand2) {
      registry.fill(HIST("hSelectionStatus"), 17);
      if ((!isDCand1) && (!isDbarCand2)) {
        registry.fill(HIST("hSelectionStatus"), 21);
      }
    }
    // Fill True info
    if (isTrueDCand1) {
      registry.fill(HIST("hSelectionStatus"), 6);
    } else if (isTrueDbarCand1) {
      registry.fill(HIST("hSelectionStatus"), 7);
    }
    if (isTrueDCand2) {
      registry.fill(HIST("hSelectionStatus"), 12);
    } else if (isTrueDbarCand2) {
      registry.fill(HIST("hSelectionStatus"), 13);
    }
    if (isTrueDCand1 && isTrueDCand2) {
      registry.fill(HIST("hSelectionStatus"), 22);
    } else if (isTrueDbarCand1 && isTrueDbarCand2) {
      registry.fill(HIST("hSelectionStatus"), 23);
    } else if (isTrueDCand1 && isTrueDbarCand2) {
      registry.fill(HIST("hSelectionStatus"), 24);
    } else if (isTrueDbarCand1 && isTrueDCand2) {
      registry.fill(HIST("hSelectionStatus"), 25);
    }
  }

  uint8_t setPairType(const uint8_t candidateType1, const uint8_t candidateType2)
  {
    uint8_t pairType(0);
    bool const isDCand1 = isD(candidateType1), isDCand2 = isD(candidateType2);
    bool const isDbarCand1 = isDbar(candidateType1), isDbarCand2 = isDbar(candidateType2);

    if (isDCand1 && isDCand2) {
      SETBIT(pairType, DD);
    }
    if (isDbarCand1 && isDbarCand2) {
      SETBIT(pairType, DbarDbar);
    }
    if (isDCand1 && isDbarCand2) {
      SETBIT(pairType, DDbar);
    }
    if (isDbarCand1 && isDCand2) {
      SETBIT(pairType, DbarD);
    }
    return pairType;
  }

  template <bool IsGen>
  void fillMatchingHisto(const int8_t matched1, const int8_t matched2)
  {
    auto hist = IsGen ? registry.get<TH1>(HIST("hMatchingMcGen")) : registry.get<TH1>(HIST("hMatchingMcRec"));
    hist->Fill(1);
    if (matched1 == DecayChannel::D0ToPiK) {
      hist->Fill(2);
    } else if (matched1 == -DecayChannel::D0ToPiK) {
      hist->Fill(3);
    } else if (matched1 == 0) {
      hist->Fill(4);
    }
    hist->Fill(5);
    if (matched2 == DecayChannel::D0ToPiK) {
      hist->Fill(6);
    } else if (matched2 == -DecayChannel::D0ToPiK) {
      hist->Fill(7);
    } else if (matched2 == 0) {
      hist->Fill(8);
    }
  }

  /// D0(bar)-D0(bar) correlation pair builder, shared by the data and MC-reco modes
  /// \param collision is the collision
  /// \param selectedD0CandidatesGrouped is the per-collision slice of the selected D0 candidates
  /// \param tracks are the (already filtered) tracks of the collision
  /// \param timeStamp is the BC timestamp of the collision
  template <bool IsMcRec, typename TCollision, typename TCandidates>
  void processImpl(TCollision const& collision,
                   TCandidates const& selectedD0CandidatesGrouped,
                   const int64_t timeStamp)
  {
    int const gCollisionId = collision.globalIndex();
    int const poolBin = corrBinning->getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));

    getCountersPerEvent<IsMcRec>(selectedD0CandidatesGrouped);

    for (const auto& candidate1 : selectedD0CandidatesGrouped) {
      analysePid(candidate1);

      outputMlD0Cand1.clear();
      outputMlD0barCand1.clear();

      auto const prong0Cand1 = candidate1.prong0Id();
      auto const prong1Cand1 = candidate1.prong1Id();

      [[maybe_unused]] int8_t matchedRec1 = 0;
      [[maybe_unused]] int8_t originRec1 = 0;
      if constexpr (IsMcRec) {
        matchedRec1 = candidate1.flagMcMatchRec();
        originRec1 = candidate1.originMcRec();
      }

      auto const candidateType1 = assignCandidateTypeD0<IsMcRec>(candidate1); // Candidate type attribution
      bool const isDCand1 = isD(candidateType1);
      bool const isDbarCand1 = isDbar(candidateType1);
      bool const isTrueDCand1 = isTrueD(candidateType1);
      bool const isTrueDbarCand1 = isTrueDbar(candidateType1);

      registry.fill(HIST("hPtCand"), candidate1.pt());

      if (!passKinematicSelections(candidate1, isDCand1, isDbarCand1)) {
        continue;
      }

      if (applyMl && !passMlSelections(candidate1, isDCand1, isDbarCand1, outputMlD0Cand1, outputMlD0barCand1)) {
        continue;
      }

      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), candidate1.y(MassD0));
      registry.fill(HIST("hPtCandAfterCut"), candidate1.pt());
      registry.fill(HIST("hPVContrib"), collision.numContrib());
      if (applyMixedEvent) {
        registry.fill(HIST("hD0Bin"), poolBin);
      }

      if (isDCand1) {
        registry.fill(HIST("hMass"), HfHelper::invMassD0ToPiK(candidate1), candidate1.pt());
        // Fill THnSparse for D0 candidates
        if (applyMl) {
          registry.fill(HIST("hnDMesonMl"), outputMlD0Cand1[0], outputMlD0Cand1[1], HfHelper::invMassD0ToPiK(candidate1), candidate1.pt(), candidate1.y(MassD0), collision.numContrib(), originRec1, candidateType1);
        } else {
          registry.fill(HIST("hnDMeson"), HfHelper::invMassD0ToPiK(candidate1), candidate1.pt(), candidate1.y(MassD0), collision.numContrib(), originRec1, candidateType1);
        }
        // Fill MC-matched histograms for D0 candidates
        if (isTrueDCand1) {
          registry.fill(HIST("hMassSig"), HfHelper::invMassD0ToPiK(candidate1), candidate1.pt());
          registry.fill(HIST("hPtVsYVsNContribMcRec"), candidate1.pt(), HfHelper::yD0(candidate1), collision.numContrib());
          registry.fill(HIST("hNContribMcRec"), collision.numContrib());
          if (originRec1 == RecoDecay::Prompt) {
            registry.fill(HIST("hMassMcRecPrompt"), HfHelper::invMassD0ToPiK(candidate1), candidate1.pt());
            registry.fill(HIST("hPtVsYVsNContribMcRecPrompt"), candidate1.pt(), HfHelper::yD0(candidate1), collision.numContrib());
          } else if (originRec1 == RecoDecay::NonPrompt) {
            registry.fill(HIST("hMassMcRecNonPrompt"), HfHelper::invMassD0ToPiK(candidate1), candidate1.pt());
            registry.fill(HIST("hPtVsYVsNContribMcRecNonPrompt"), candidate1.pt(), HfHelper::yD0(candidate1), collision.numContrib());
          }
        } else if (isTrueDbarCand1) {
          registry.fill(HIST("hMassMcRecReflections"), HfHelper::invMassD0ToPiK(candidate1), candidate1.pt());
        }
        // Fill D0 table for offline event mixing
        if (applyMixedEvent) {
          entryDMesonInfoCand1(candidate1.pt(), candidate1.eta(), candidate1.phi(), HfHelper::invMassD0ToPiK(candidate1), candidateType1, poolBin, gCollisionId, timeStamp);
        }
      }
      if (isDbarCand1) {
        registry.fill(HIST("hMass"), HfHelper::invMassD0barToKPi(candidate1), candidate1.pt());
        // Fill THnSparse for D0bar candidates
        if (applyMl) {
          registry.fill(HIST("hnDMesonMl"), outputMlD0barCand1[0], outputMlD0barCand1[1], HfHelper::invMassD0barToKPi(candidate1), candidate1.pt(), candidate1.y(MassD0), collision.numContrib(), originRec1, candidateType1);
        } else {
          registry.fill(HIST("hnDMeson"), HfHelper::invMassD0barToKPi(candidate1), candidate1.pt(), candidate1.y(MassD0), collision.numContrib(), originRec1, candidateType1);
        }
        // Fill MC-matched histograms for D0bar candidates
        if (isTrueDbarCand1) {
          registry.fill(HIST("hMassSig"), HfHelper::invMassD0barToKPi(candidate1), candidate1.pt());
          registry.fill(HIST("hPtVsYVsNContribMcRec"), candidate1.pt(), HfHelper::yD0(candidate1), collision.numContrib());
          registry.fill(HIST("hNContribMcRec"), collision.numContrib());
          if (originRec1 == RecoDecay::Prompt) {
            registry.fill(HIST("hMassMcRecPrompt"), HfHelper::invMassD0barToKPi(candidate1), candidate1.pt());
            registry.fill(HIST("hPtVsYVsNContribMcRecPrompt"), candidate1.pt(), HfHelper::yD0(candidate1), collision.numContrib());
          } else if (originRec1 == RecoDecay::NonPrompt) {
            registry.fill(HIST("hMassMcRecNonPrompt"), HfHelper::invMassD0barToKPi(candidate1), candidate1.pt());
            registry.fill(HIST("hPtVsYVsNContribMcRecNonPrompt"), candidate1.pt(), HfHelper::yD0(candidate1), collision.numContrib());
          }
        } else if (isTrueDCand1) {
          registry.fill(HIST("hMassMcRecReflections"), HfHelper::invMassD0barToKPi(candidate1), candidate1.pt());
        }
        // Fill D0bar table for offline event mixing
        if (applyMixedEvent) {
          entryDMesonInfoCand1(candidate1.pt(), candidate1.eta(), candidate1.phi(), HfHelper::invMassD0barToKPi(candidate1), candidateType1, poolBin, gCollisionId, timeStamp);
        }
      }

      // Loop on the second D0 candidate
      for (auto candidate2 = candidate1 + 1; candidate2 != selectedD0CandidatesGrouped.end(); ++candidate2) {

        outputMlD0Cand2.clear();
        outputMlD0barCand2.clear();

        auto const prong0Cand2 = candidate2.prong0Id();
        auto const prong1Cand2 = candidate2.prong1Id();
        if ((prong0Cand1 == prong0Cand2) || (prong1Cand1 == prong1Cand2) || (prong0Cand1 == prong1Cand2) || (prong1Cand1 == prong0Cand2)) {
          continue;
        }

        auto const candidateType2 = assignCandidateTypeD0<IsMcRec>(candidate2); // Candidate type attribution
        bool const isDCand2 = isD(candidateType2);
        bool const isDbarCand2 = isDbar(candidateType2);

        int8_t matchedRec2 = 0;
        int8_t originRec2 = 0;
        if constexpr (IsMcRec) {
          matchedRec2 = candidate2.flagMcMatchRec();
          originRec2 = candidate2.originMcRec();
        }

        if (!passKinematicSelections(candidate2, isDCand2, isDbarCand2)) {
          continue;
        }

        if (applyMl && !passMlSelections(candidate2, isDCand2, isDbarCand2, outputMlD0Cand2, outputMlD0barCand2)) {
          continue;
        }

        if (applyMixedEvent) {
          if (isDCand2) {
            entryDMesonInfoCand2(candidate2.pt(), candidate2.eta(), candidate2.phi(), HfHelper::invMassD0ToPiK(candidate2), candidateType2, poolBin, gCollisionId, timeStamp);
          }
          if (isDbarCand2) {
            entryDMesonInfoCand2(candidate2.pt(), candidate2.eta(), candidate2.phi(), HfHelper::invMassD0barToKPi(candidate2), candidateType2, poolBin, gCollisionId, timeStamp);
          }
          registry.fill(HIST("hD0Cand2Bin"), poolBin);
        }

        // Canonical ordering of the pair: D0 first (higher pT first for equal flavours)
        bool const swapOrder = shouldSwapOrder(isDCand1, isDbarCand1, isDCand2, isDbarCand2, candidate1, candidate2);
        uint8_t const type1 = swapOrder ? candidateType2 : candidateType1;
        uint8_t const type2 = swapOrder ? candidateType1 : candidateType2;
        auto const& first = swapOrder ? candidate2 : candidate1;
        auto const& second = swapOrder ? candidate1 : candidate2;
        uint8_t const pairType = setPairType(type1, type2);

        // Fill entries for D0-D0(bar) pairs and plot histos
        fillSelectionStatusHisto(type1, type2);
        entryD0Pair(first.pt(), second.pt(), HfHelper::yD0(first), HfHelper::yD0(second), first.eta(), second.eta(),
                    first.phi(), second.phi(), HfHelper::invMassD0ToPiK(first), HfHelper::invMassD0barToKPi(first),
                    HfHelper::invMassD0ToPiK(second), HfHelper::invMassD0barToKPi(second), pairType, type1, type2);
        if (applyMl) {
          entryD0PairMl(swapOrder ? outputMlD0Cand2 : outputMlD0Cand1, swapOrder ? outputMlD0barCand2 : outputMlD0barCand1,
                        swapOrder ? outputMlD0Cand1 : outputMlD0Cand2, swapOrder ? outputMlD0barCand1 : outputMlD0barCand2);
        }
        if constexpr (IsMcRec) {
          fillMatchingHisto<false>(swapOrder ? matchedRec2 : matchedRec1, swapOrder ? matchedRec1 : matchedRec2);
          entryD0PairMcInfo(swapOrder ? originRec2 : originRec1, swapOrder ? originRec1 : originRec2,
                            swapOrder ? matchedRec2 : matchedRec1, swapOrder ? matchedRec1 : matchedRec2);
        }
      } // end inner loop (Cand2)
    } // end outer loop (Cand1)

    if (applyMixedEvent) {
      registry.fill(HIST("hZvtx"), collision.posZ());
      registry.fill(HIST("hMultFT0M"), collision.multFT0M());
    }
  }

  /// D0(bar)-D0(bar) correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(SelCollisionsWithD0::iterator const& collision,
                   CandidatesD0Data const&, aod::Tracks const&, aod::BCsWithTimestamps const&)
  {
    auto const timeStamp = collision.bc_as<aod::BCsWithTimestamps>().timestamp();
    auto grouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    processImpl<false>(collision, grouped, timeStamp);
  }
  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processData, "Process data mode", true);

  /// D0(bar)-D0(bar) correlation pair builder - MC reco level (candidates matched to MC truth)
  void processMcRec(SelCollisionsWithD0::iterator const& collision,
                    CandidatesD0McRec const&, aod::Tracks const&, aod::BCsWithTimestamps const&)
  {
    auto const timeStamp = collision.bc_as<aod::BCsWithTimestamps>().timestamp();
    auto grouped = selectedD0CandidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    processImpl<true>(collision, grouped, timeStamp);
  }
  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processMcRec, "Process MC reco mode", false);

  /// D0(bar)-D0(bar) correlation pair builder - MC generator level (particles flagged as D0 -> pi K)
    void processMcGen(aod::McCollision const&, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, McParticlesPlus2Prong const& mcParticles)
  {
    // Take the associated reconstructed collision with the highest number of PV contributors
    int numPvContributorsGen{0};
    for (const auto& collision : collisions) {
      int const numPvContributors = collision.numContrib();
      if (numPvContributors > numPvContributorsGen) {
        numPvContributorsGen = numPvContributors;
      }
    }

    // Select the D0(bar) -> K pi particles
    std::vector<McParticlesPlus2Prong::iterator> selectedParticles;
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) != DecayChannel::D0ToPiK) {
        continue;
      }
      registry.fill(HIST("hPtCandMcGen"), particle.pt());
      if (!passGenSelections(particle)) {
        continue;
      }
      selectedParticles.push_back(particle);
    }

    // Single-particle histograms and per-event counters
    int nDevent = 0, nDbarevent = 0, nDDbarevent = 0, nDorDbarevent = 0;
    for (const auto& particle : selectedParticles) {
      auto const particleType = assignParticleTypeD0Gen(particle);
      bool const isDParticle = isTrueD(particleType);
      bool const isDbarParticle = isTrueDbar(particleType);
      nDevent += isDParticle;
      nDbarevent += isDbarParticle;
      nDDbarevent += (isDParticle && isDbarParticle);
      nDorDbarevent += (isDParticle || isDbarParticle);

      registry.fill(HIST("hEtaMcGen"), particle.eta());
      registry.fill(HIST("hPhiMcGen"), particle.phi());
      registry.fill(HIST("hPtCandAfterCutMcGen"), particle.pt());

      registry.fill(HIST("hStatusSinglePartMcGen"), 1);
      if (isDParticle && !isDbarParticle) {
        registry.fill(HIST("hStatusSinglePartMcGen"), 2);
        registry.fill(HIST("hStatusSinglePartMcGen"), 5);
      } else if (isDbarParticle && !isDParticle) {
        registry.fill(HIST("hStatusSinglePartMcGen"), 3);
        registry.fill(HIST("hStatusSinglePartMcGen"), 6);
      } else if (isDParticle && isDbarParticle) {
        registry.fill(HIST("hStatusSinglePartMcGen"), 4);
      }

      registry.fill(HIST("hPtVsYVsNContribMcGen"), particle.pt(), particle.y(), numPvContributorsGen);
      if (particle.originMcGen() == RecoDecay::Prompt) {
        registry.fill(HIST("hPtVsYVsNContribMcGenPrompt"), particle.pt(), particle.y(), numPvContributorsGen);
      } else if (particle.originMcGen() == RecoDecay::NonPrompt) {
        registry.fill(HIST("hPtVsYVsNContribMcGenNonPrompt"), particle.pt(), particle.y(), numPvContributorsGen);
      }
      registry.fill(HIST("hNContribMcGen"), numPvContributorsGen);
    }
    if (nDevent > 0) {
      registry.fill(HIST("hInputCheckD0McGen"), nDevent);
    }
    if (nDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0barMcGen"), nDbarevent);
    }
    if (nDDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0AndD0barMcGen"), nDDbarevent);
    }
    if (nDorDbarevent > 0) {
      registry.fill(HIST("hInputCheckD0OrD0barMcGen"), nDorDbarevent);
    }

    // Pairs of selected particles
    for (size_t i = 0; i < selectedParticles.size(); ++i) {
      for (size_t j = i + 1; j < selectedParticles.size(); ++j) {
        auto const& particle1 = selectedParticles[i];
        auto const& particle2 = selectedParticles[j];

        auto const type1 = assignParticleTypeD0Gen(particle1);
        auto const type2 = assignParticleTypeD0Gen(particle2);

        // Canonical ordering (same as reco): D0 first, higher pT first for equal flavours
        bool const swapOrder = shouldSwapOrder(isTrueD(type1), isTrueDbar(type1), isTrueD(type2), isTrueDbar(type2), particle1, particle2);
        auto const& first = swapOrder ? particle2 : particle1;
        auto const& second = swapOrder ? particle1 : particle2;
        auto const typeFirst = swapOrder ? type2 : type1;
        auto const typeSecond = swapOrder ? type1 : type2;
        bool const isDFirst = isTrueD(typeFirst), isDbarFirst = isTrueDbar(typeFirst);
        bool const isDSecond = isTrueD(typeSecond), isDbarSecond = isTrueDbar(typeSecond);

        fillMatchingHisto<true>(first.flagMcMatchGen(), second.flagMcMatchGen());

        // Selection status
        registry.fill(HIST("hSelectionStatusMcGen"), 2);
        if (isDFirst && !isDbarFirst) {
          registry.fill(HIST("hSelectionStatusMcGen"), 6);
        } else if (isDbarFirst && !isDFirst) {
          registry.fill(HIST("hSelectionStatusMcGen"), 7);
        }
        registry.fill(HIST("hSelectionStatusMcGen"), 8);
        if (isDSecond && !isDbarSecond) {
          registry.fill(HIST("hSelectionStatusMcGen"), 12);
        } else if (isDbarSecond && !isDSecond) {
          registry.fill(HIST("hSelectionStatusMcGen"), 13);
        }

        registry.fill(HIST("hSelectionStatusMcGen"), 1);
        if (isDFirst && isDSecond) {
          registry.fill(HIST("hSelectionStatusMcGen"), 22);
        }
        if (isDbarFirst && isDbarSecond) {
          registry.fill(HIST("hSelectionStatusMcGen"), 23);
        }
        if (isDFirst && isDbarSecond) {
          registry.fill(HIST("hSelectionStatusMcGen"), 24);
        }
        if (isDbarFirst && isDSecond) {
          registry.fill(HIST("hSelectionStatusMcGen"), 25);
        }

        uint8_t const pairType = setPairType(typeFirst, typeSecond);
        entryD0PairMcGen(first.pt(), second.pt(), first.y(), second.y(), first.eta(), second.eta(), first.phi(), second.phi(),
                         MassD0, MassD0Bar, MassD0, MassD0Bar, pairType, typeFirst, typeSecond);
        entryD0PairMcGenInfo(first.originMcGen(), second.originMcGen(), first.flagMcMatchGen(), second.flagMcMatchGen());
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processMcGen, "Process D0 Mc Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCorrelatorD0HadronsSelection>(cfgc),
    adaptAnalysisTask<HfCorrelatorDMesonPairs>(cfgc),
  };
}
