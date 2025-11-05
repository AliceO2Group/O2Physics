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

/// \file taskD0.cxx
/// \brief D0 analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>

#include <algorithm> // std::min
#include <array>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;

/// D0 analysis task
namespace
{
enum CandTypeSel {
  SigD0 = 0,      // Signal D0
  SigD0bar,       // Signal D0bar
  ReflectedD0,    // Reflected D0
  ReflectedD0bar, // Reflected D0bar
  PureSigD0,      // Signal D0 exclude Reflected D0bar
  PureSigD0bar    // Signal D0bar exclude Reflected D0
};
} // namespace
struct HfTaskD0 {
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<int> selectionFlagHf{"selectionFlagHf", 1, "Selection Flag for HF flagged candidates"};
  Configurable<int> selectionTopol{"selectionTopol", 1, "Selection Flag for topologically selected candidates"};
  Configurable<int> selectionCand{"selectionCand", 1, "Selection Flag for conj. topol. selected candidates"};
  Configurable<int> selectionPid{"selectionPid", 1, "Selection Flag for reco PID candidates"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  Configurable<int> centEstimator{"centEstimator", 0, "Centrality estimation (None: 0, FT0C: 2, FT0M: 3)"};
  Configurable<int> occEstimator{"occEstimator", 0, "Occupancy estimation (None: 0, ITS: 1, FT0C: 2)"};
  Configurable<bool> storeCentrality{"storeCentrality", false, "Flag to store centrality information"};
  Configurable<bool> storeOccupancyAndIR{"storeOccupancyAndIR", false, "Flag to store occupancy information and interaction rate"};
  Configurable<bool> storeTrackQuality{"storeTrackQuality", false, "Flag to store track quality information"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> irSource{"irSource", "ZNC hadronic", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  ctpRateFetcher mRateFetcher;

  SliceCache cache;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using D0Candidates = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
  using D0CandidatesMc = soa::Join<D0Candidates, aod::HfCand2ProngMcRec>;
  using D0CandidatesKF = soa::Join<D0Candidates, aod::HfCand2ProngKF>;
  using D0CandidatesMcKF = soa::Join<D0CandidatesKF, aod::HfCand2ProngMcRec>;

  using D0CandidatesMl = soa::Join<D0Candidates, aod::HfMlD0>;
  using D0CandidatesMlMc = soa::Join<D0CandidatesMl, aod::HfCand2ProngMcRec>;
  using D0CandidatesMlKF = soa::Join<D0CandidatesMl, aod::HfCand2ProngKF>;
  using D0CandidatesMlMcKF = soa::Join<D0CandidatesMlKF, aod::HfCand2ProngMcRec>;

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs>;
  using CollisionsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using CollisionsWithMcLabelsCent = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs>;
  using TracksSelQuality = soa::Join<aod::TracksExtra, aod::TracksWMc>;
  PresliceUnsorted<CollisionsWithMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<CollisionsWithMcLabelsCent> colPerMcCollisionCent = aod::mccollisionlabel::mcCollisionId;

  Partition<D0Candidates> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesKF> selectedD0CandidatesKF = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesMc> selectedD0CandidatesMc = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;
  Partition<D0CandidatesMcKF> selectedD0CandidatesMcKF = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;

  Partition<D0CandidatesMl> selectedD0CandidatesMl = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesMlKF> selectedD0CandidatesMlKF = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesMlMc> selectedD0CandidatesMlMc = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;
  Partition<D0CandidatesMlMcKF> selectedD0CandidatesMlMcKF = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;

  // ThnSparse for ML outputScores and Vars
  ConfigurableAxis thnConfigAxisBkgScore{"thnConfigAxisBkgScore", {50, 0, 1}, "Bkg score bins"};
  ConfigurableAxis thnConfigAxisNonPromptScore{"thnConfigAxisNonPromptScore", {50, 0, 1}, "Non-prompt score bins"};
  ConfigurableAxis thnConfigAxisPromptScore{"thnConfigAxisPromptScore", {50, 0, 1}, "Prompt score bins"};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {120, 1.5848, 2.1848}, "Cand. inv-mass bins"};
  ConfigurableAxis thnConfigAxisPtB{"thnConfigAxisPtB", {1000, 0, 100}, "Cand. beauty mother pTB bins"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {500, 0, 50}, "Cand. pT bins"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, "Cand. origin type"};
  ConfigurableAxis thnConfigAxisCandType{"thnConfigAxisCandType", {6, -0.5, 5.5}, "D0 type"};
  ConfigurableAxis thnConfigAxisGenPtD{"thnConfigAxisGenPtD", {500, 0, 50}, "Gen Pt D"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {1000, 0, 100}, "Gen Pt B"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {110, 0., 110.}, ""};
  ConfigurableAxis thnConfigAxisOccupancy{"thnConfigAxisOccupancy", {14, 0, 14000}, "axis for centrality"};
  ConfigurableAxis thnConfigAxisMinItsNCls{"thnConfigAxisMinItsNCls", {5, 3, 8}, "axis for minimum ITS NCls of candidate prongs"};
  ConfigurableAxis thnConfigAxisMinTpcNCrossedRows{"thnConfigAxisMinTpcNCrossedRows", {10, 70, 180}, "axis for minimum TPC NCls crossed rows of candidate prongs"};
  ConfigurableAxis thnConfigAxisIR{"thnConfigAxisIR", {5000, 0, 500}, "Interaction rate (kHz)"};

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSig", "2-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSigPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSigNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecBg", "2-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hYGenPrompt", "MC particles (matched, prompt);#it{y}^{gen.};entries", {HistType::kTH1F, {{300, -1.5, 1.5}}}},
     {"hYGenNonPrompt", "MC particles (matched, non-prompt);#it{y}^{gen.};entries", {HistType::kTH1F, {{300, -1.5, 1.5}}}},
     {"hPtGenSig", "2-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hCPARecSig", "2-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecBg", "2-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hEtaRecSig", "2-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hEtaRecBg", "2-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hPtGenVsPtRecSig", "2-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH2F, {{360, 0., 36.}, {360, 0., 36.}}}},
     {"hYGenVsYRecSig", "2-prong candidates (matched);#it{y}^{gen.} ;#it{y}^{rec.} ;entries", {HistType::kTH2F, {{300, -1.5, 1.5}, {300, -1.5, 1.5}}}},
     {"hPtProng0Sig", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hPtProng1Sig", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hDecLengthSig", "2-prong candidates (matched);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hDecLengthXYSig", "2-prong candidates (matched);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthSig", "2-prong candidates (matched);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthXYSig", "2-prong candidates (matched);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hd0Prong0Sig", "2-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0Prong1Sig", "2-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0d0Sig", "2-prong candidates (matched);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}}},
     {"hCTSSig", "2-prong candidates (matched);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCtSig", "2-prong candidates (matched);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}}},
     {"hCPASig", "2-prong candidates (matched);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCPAxySig", "2-prong candidates (matched);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hPtProng0Bkg", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hPtProng1Bkg", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hDecLengthBkg", "2-prong candidates (checked);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hDecLengthXYBkg", "2-prong candidates (checked);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthBkg", "2-prong candidates (checked);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthXYBkg", "2-prong candidates (checked);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hd0Prong0Bkg", "2-prong candidates (checked);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0Prong1Bkg", "2-prong candidates (checked);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0d0Bkg", "2-prong candidates (checked);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}}},
     {"hCTSBkg", "2-prong candidates (checked);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCtBkg", "2-prong candidates (checked);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}}},
     {"hCPABkg", "2-prong candidates (checked);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCPAxyBkg", "2-prong candidates (checked);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hPtVsYRecSig_RecoPID", "2-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoPID", "2-prong candidates (RecoPID - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoPID", "2-prong candidates (RecoPID - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoCand", "2-prong candidates (RecoCand - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoCand", "2-prong candidates (RecoCand - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoCand", "2-prong candidates (RecoCand - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoTopol", "2-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoTopol", "2-prong candidates (RecoTopol - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoTopol", "2-prong candidates (RecoTopol - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoHFFlag", "2-prong candidates (RecoHFFlag - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoHFFlag", "2-prong candidates (RecoHFFlag - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoHFFlag", "2-prong candidates (RecoHFFlag - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGen", "2-prong candidates (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGenPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGenNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hMassVsPtGenVsPtRecSig", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{gen.} (GeV/#it{c});#it{p}_{T}^{rec.} (GeV/#it{c})", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {150, 0., 30.}}}},
     {"hMassSigD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassBkgD0", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassReflBkgD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigBkgD0", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassBkgD0bar", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassReflBkgD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigBkgD0bar", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}}}};

  void init(InitContext&)
  {
    std::array<bool, 12> doprocess{doprocessDataWithDCAFitterN, doprocessDataWithDCAFitterNCent, doprocessDataWithKFParticle, doprocessMcWithDCAFitterN, doprocessMcWithDCAFitterNCent, doprocessMcWithKFParticle, doprocessDataWithDCAFitterNMl, doprocessDataWithDCAFitterNMlCent, doprocessDataWithKFParticleMl, doprocessMcWithDCAFitterNMl, doprocessMcWithDCAFitterNMlCent, doprocessMcWithKFParticleMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) == 0) {
      LOGP(fatal, "At least one process function should be enabled at a time.");
    }
    if ((doprocessDataWithDCAFitterN || doprocessDataWithDCAFitterNCent || doprocessMcWithDCAFitterN || doprocessMcWithDCAFitterNCent || doprocessDataWithDCAFitterNMl || doprocessDataWithDCAFitterNMlCent || doprocessMcWithDCAFitterNMl || doprocessMcWithDCAFitterNMlCent) && (doprocessDataWithKFParticle || doprocessMcWithKFParticle || doprocessDataWithKFParticleMl || doprocessMcWithKFParticleMl)) {
      LOGP(fatal, "DCAFitterN and KFParticle can not be enabled at a time.");
    }
    if ((storeCentrality || storeOccupancyAndIR) && !(doprocessDataWithDCAFitterNCent || doprocessMcWithDCAFitterNCent || doprocessDataWithDCAFitterNMlCent || doprocessMcWithDCAFitterNMlCent)) {
      LOGP(fatal, "Can't enable the storeCentrality and storeOccupancu without cent process");
    }
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassVsPhi", "2-prong candidates vs phi;inv. mass (#pi K) (GeV/#it{c}^{2});phi (rad);entries", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {32, 0, o2::constants::math::TwoPI}}});
    registry.add("hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0ErrProng0", "2-prong candidates;prong 0 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0ErrProng1", "2-prong candidates;prong 1 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassFinerBinning", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthFinerBinning", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{400, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxyFinerBinning", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{400, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0FinerBinning", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1FinerBinning", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0FinerBinning", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -0.1, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTSFinerBinning", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtFinerBinning", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{500, -0., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAFinerBinning", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAXYFinerBinning", "2-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0VsPtSig", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1VsPtSig", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0VsPtSig", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2}) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -0.1, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTSVsPtSig", "2-prong candidates;cos #it{#theta}* (D^{0}) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAVsPtSig", "2-prong candidates;cosine of pointing angle vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAXYVsPtSig", "2-prong candidates;cosine of pointing angle xy vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLengthxyVsPtSig", "2-prong candidates;decay length xy (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthVsPtSig", "2-prong candidates;decay length (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxyVsPtSig", "2-prong candidates;decay length xy (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (#pi K) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPtB{thnConfigAxisPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisY{thnConfigAxisY, "y"};
    const AxisSpec thnAxisOrigin{thnConfigAxisOrigin, "Origin"};
    const AxisSpec thnAxisCandType{thnConfigAxisCandType, "D0 type"};
    const AxisSpec thnAxisGenPtD{thnConfigAxisGenPtD, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisGenPtB{thnConfigAxisGenPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisNumPvContr{thnConfigAxisNumPvContr, "Number of PV contributors"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality"};
    const AxisSpec thnAxisOccupancy{thnConfigAxisOccupancy, "Occupancy"};
    const AxisSpec thnAxisMinItsNCls{thnConfigAxisMinItsNCls, "Minimum ITS cluster found"};
    const AxisSpec thnAxisMinTpcNCrossedRows{thnConfigAxisMinTpcNCrossedRows, "Minimum TPC crossed rows"};
    const AxisSpec thnAxisIR{thnConfigAxisIR, "Interaction rate"};

    if (doprocessMcWithDCAFitterN || doprocessMcWithDCAFitterNCent || doprocessMcWithKFParticle || doprocessMcWithDCAFitterNMl || doprocessMcWithDCAFitterNMlCent || doprocessMcWithKFParticleMl) {
      std::vector<AxisSpec> axesAcc = {thnAxisGenPtD, thnAxisGenPtB, thnAxisY, thnAxisOrigin, thnAxisNumPvContr};

      if (storeCentrality) {
        axesAcc.push_back(thnAxisCent);
      }
      // interaction rate only store in Data and MC Reco. Level
      if (storeOccupancyAndIR) {
        axesAcc.push_back(thnAxisOccupancy);
      }

      registry.add("hSparseAcc", "Thn for generated D0 from charm and beauty", HistType::kTHnSparseD, axesAcc);
      registry.get<THnSparse>(HIST("hSparseAcc"))->Sumw2();
    }

    std::vector<AxisSpec> axes = {thnAxisMass, thnAxisPt, thnAxisY, thnAxisCandType};
    if (doprocessMcWithDCAFitterN || doprocessMcWithDCAFitterNCent || doprocessMcWithKFParticle || doprocessMcWithDCAFitterNMl || doprocessMcWithDCAFitterNMlCent || doprocessMcWithKFParticleMl) {
      axes.push_back(thnAxisPtB);
      axes.push_back(thnAxisOrigin);
      axes.push_back(thnAxisNumPvContr);
    }
    if (storeCentrality) {
      axes.push_back(thnAxisCent);
    }
    if (storeOccupancyAndIR) {
      axes.push_back(thnAxisOccupancy);
      axes.push_back(thnAxisIR);
    }
    if (storeTrackQuality) {
      axes.push_back(thnAxisMinItsNCls);
      axes.push_back(thnAxisMinTpcNCrossedRows);
    }
    if (applyMl) {
      const AxisSpec thnAxisBkgScore{thnConfigAxisBkgScore, "BDT score bkg."};
      const AxisSpec thnAxisNonPromptScore{thnConfigAxisNonPromptScore, "BDT score non-prompt."};
      const AxisSpec thnAxisPromptScore{thnConfigAxisPromptScore, "BDT score prompt."};

      axes.insert(axes.begin(), thnAxisPromptScore);
      axes.insert(axes.begin(), thnAxisNonPromptScore);
      axes.insert(axes.begin(), thnAxisBkgScore);

      registry.add("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type", "Thn for D0 candidates", HistType::kTHnSparseD, axes);
      registry.get<THnSparse>(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"))->Sumw2();
    } else {
      registry.add("hMassVsPtVsPtBVsYVsOriginVsD0Type", "Thn for D0 candidates", HistType::kTHnSparseD, axes);
      registry.get<THnSparse>(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"))->Sumw2();
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  template <int ReconstructionType, bool ApplyMl, typename CandType, typename CollType, typename BCsType>
  void processData(CandType const& candidates,
                   CollType const&,
                   aod::TracksWExtra const&,
                   BCsType const&)
  {
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yD0(candidate)) > yCandRecoMax) {
        continue;
      }

      float massD0, massD0bar;
      if constexpr (ReconstructionType == aod::hf_cand::VertexerType::KfParticle) {
        massD0 = candidate.kfGeoMassD0();
        massD0bar = candidate.kfGeoMassD0bar();
      } else {
        massD0 = HfHelper::invMassD0ToPiK(candidate);
        massD0bar = HfHelper::invMassD0barToKPi(candidate);
      }
      auto ptCandidate = candidate.pt();

      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), massD0, ptCandidate);
        registry.fill(HIST("hMassFinerBinning"), massD0, ptCandidate);
        registry.fill(HIST("hMassVsPhi"), massD0, ptCandidate, candidate.phi());
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), massD0bar, ptCandidate);
        registry.fill(HIST("hMassFinerBinning"), massD0bar, ptCandidate);
        registry.fill(HIST("hMassVsPhi"), massD0bar, ptCandidate, candidate.phi());
      }
      registry.fill(HIST("hPtCand"), ptCandidate);
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), ptCandidate);
      registry.fill(HIST("hDecLengthxy"), candidate.decayLengthXY(), ptCandidate);
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), ptCandidate);
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), ptCandidate);
      registry.fill(HIST("hNormalisedDecLength"), candidate.decayLengthNormalised(), ptCandidate);
      registry.fill(HIST("hNormalisedDecLengthxy"), candidate.decayLengthXYNormalised(), ptCandidate);
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), ptCandidate);
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), ptCandidate);
      registry.fill(HIST("hd0ErrProng0"), candidate.errorImpactParameter0(), ptCandidate);
      registry.fill(HIST("hd0ErrProng1"), candidate.errorImpactParameter1(), ptCandidate);
      registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), ptCandidate);
      registry.fill(HIST("hCTS"), HfHelper::cosThetaStarD0(candidate), ptCandidate);
      registry.fill(HIST("hCt"), HfHelper::ctD0(candidate), ptCandidate);
      registry.fill(HIST("hCPA"), candidate.cpa(), ptCandidate);
      registry.fill(HIST("hEta"), candidate.eta(), ptCandidate);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), ptCandidate);
      registry.fill(HIST("hDecLengthFinerBinning"), candidate.decayLength(), ptCandidate);
      registry.fill(HIST("hDecLengthxyFinerBinning"), candidate.decayLengthXY(), ptCandidate);
      registry.fill(HIST("hd0Prong0FinerBinning"), candidate.impactParameter0(), ptCandidate);
      registry.fill(HIST("hd0Prong1FinerBinning"), candidate.impactParameter1(), ptCandidate);
      registry.fill(HIST("hd0d0FinerBinning"), candidate.impactParameterProduct(), ptCandidate);
      registry.fill(HIST("hCTSFinerBinning"), HfHelper::cosThetaStarD0(candidate), ptCandidate);
      registry.fill(HIST("hCtFinerBinning"), HfHelper::ctD0(candidate), ptCandidate);
      registry.fill(HIST("hCPAFinerBinning"), candidate.cpa(), ptCandidate);
      registry.fill(HIST("hCPAXYFinerBinning"), candidate.cpaXY(), ptCandidate);

      float cent{-1.f};
      float occ{-1.f};
      float ir{-1.f};
      if (storeCentrality || storeOccupancyAndIR) {
        auto collision = candidate.template collision_as<CollType>();
        if (storeCentrality && centEstimator != CentralityEstimator::None) {
          cent = getCentralityColl(collision, centEstimator);
        }
        if (storeOccupancyAndIR && occEstimator != OccupancyEstimator::None) {
          occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
          auto bc = collision.template foundBC_as<BCsType>();
          ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource, true) * 1.e-3; // kHz
        }
      }

      auto trackPos = candidate.template prong0_as<o2::aod::TracksWExtra>(); // positive daughter
      auto trackNeg = candidate.template prong1_as<o2::aod::TracksWExtra>(); // negative daughter
      int const minItsClustersOfProngs = std::min(trackPos.itsNCls(), trackNeg.itsNCls());
      int const minTpcCrossedRowsOfProngs = std::min(trackPos.tpcNClsCrossedRows(), trackNeg.tpcNClsCrossedRows());
      if constexpr (ApplyMl) {
        if (storeCentrality && storeOccupancyAndIR) {
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, HfHelper::yD0(candidate), SigD0, cent, occ, ir);
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0, cent, occ, ir);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, HfHelper::yD0(candidate), SigD0bar, cent, occ, ir);
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar, cent, occ, ir);
          }
        } else if (storeCentrality && !storeOccupancyAndIR) {
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, HfHelper::yD0(candidate), SigD0, cent);
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0, cent);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, HfHelper::yD0(candidate), SigD0bar, cent);
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar, cent);
          }
        } else if (!storeCentrality && storeOccupancyAndIR) {
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, HfHelper::yD0(candidate), SigD0, occ, ir);
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0, occ, ir);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, HfHelper::yD0(candidate), SigD0bar, occ, ir);
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar, occ, ir);
          }
        } else if (storeTrackQuality) {
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, HfHelper::yD0(candidate), SigD0, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, HfHelper::yD0(candidate), SigD0bar, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
          }
        } else {
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, HfHelper::yD0(candidate), SigD0);
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, HfHelper::yD0(candidate), SigD0bar);
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar);
          }
        }
      } else {
        if (storeCentrality && storeOccupancyAndIR) {
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, HfHelper::yD0(candidate), SigD0, cent, occ, ir);
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0, cent, occ, ir);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, HfHelper::yD0(candidate), SigD0bar, cent, occ, ir);
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar, cent, occ, ir);
          }
        } else if (storeCentrality && !storeOccupancyAndIR) {
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, HfHelper::yD0(candidate), SigD0, cent);
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0, cent);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, HfHelper::yD0(candidate), SigD0bar, cent);
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar, cent);
          }
        } else if (!storeCentrality && storeOccupancyAndIR) {
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, HfHelper::yD0(candidate), SigD0, occ, ir);
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0, occ, ir);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, HfHelper::yD0(candidate), SigD0bar, occ, ir);
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar, occ, ir);
          }
        } else if (storeTrackQuality) {
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, HfHelper::yD0(candidate), SigD0, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, HfHelper::yD0(candidate), SigD0bar);
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar);
          }
        } else {
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, HfHelper::yD0(candidate), SigD0);
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, HfHelper::yD0(candidate), SigD0bar);
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar);
          }
        }
      }
    }
  }
  void processDataWithDCAFitterN(D0Candidates const&, Collisions const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter, false>(selectedD0Candidates, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processDataWithDCAFitterN, "process taskD0 with DCAFitterN", true);

  void processDataWithDCAFitterNCent(D0Candidates const&, CollisionsCent const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter, false>(selectedD0Candidates, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processDataWithDCAFitterNCent, "process taskD0 with DCAFitterN and centrality", false);

  void processDataWithKFParticle(D0CandidatesKF const&, Collisions const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::KfParticle, false>(selectedD0CandidatesKF, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processDataWithKFParticle, "process taskD0 with KFParticle", false);
  // TODO: add processKFParticleCent

  void processDataWithDCAFitterNMl(D0CandidatesMl const&, Collisions const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter, true>(selectedD0CandidatesMl, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processDataWithDCAFitterNMl, "process taskD0 with DCAFitterN and ML selections", false);

  void processDataWithDCAFitterNMlCent(D0CandidatesMl const&, CollisionsCent const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter, true>(selectedD0CandidatesMl, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processDataWithDCAFitterNMlCent, "process taskD0 with DCAFitterN and ML selections and centrality", false);

  void processDataWithKFParticleMl(D0CandidatesMlKF const&, Collisions const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::KfParticle, true>(selectedD0CandidatesMlKF, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processDataWithKFParticleMl, "process taskD0 with KFParticle and ML selections", false);
  // TODO: add processKFParticleMlCent

  template <int ReconstructionType, bool ApplyMl, typename CandType, typename CollType, typename BCsType>
  void processMc(CandType const& candidates,
                 soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
                 TracksSelQuality const&,
                 CollType const& collisions,
                 aod::McCollisions const&,
                 BCsType const&)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yD0(candidate)) > yCandRecoMax) {
        continue;
      }

      float cent{-1.f};
      float occ{-1.f};
      float ir{-1.f};
      auto collision = candidate.template collision_as<CollType>();
      auto numPvContributors = collision.numContrib();
      if (storeCentrality && centEstimator != CentralityEstimator::None) {
        cent = getCentralityColl(collision, centEstimator);
      }
      if (storeOccupancyAndIR && occEstimator != OccupancyEstimator::None) {
        occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
        auto bc = collision.template foundBC_as<BCsType>();
        ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource, true) * 1.e-3; // kHz
      }
      float massD0, massD0bar;
      if constexpr (ReconstructionType == aod::hf_cand::VertexerType::KfParticle) {
        massD0 = candidate.kfGeoMassD0();
        massD0bar = candidate.kfGeoMassD0bar();
      } else {
        massD0 = HfHelper::invMassD0ToPiK(candidate);
        massD0bar = HfHelper::invMassD0barToKPi(candidate);
      }
      auto trackPos = candidate.template prong0_as<TracksSelQuality>(); // positive daughter
      auto trackNeg = candidate.template prong1_as<TracksSelQuality>(); // negative daughter
      if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(mcParticles, trackPos.template mcParticle_as<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>(), o2::constants::physics::Pdg::kD0, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);
        auto ptGen = particleMother.pt();                                                   // gen. level pT
        auto yGen = RecoDecay::y(particleMother.pVector(), o2::constants::physics::MassD0); // gen. level y
        registry.fill(HIST("hPtGenSig"), ptGen);                                            // gen. level pT
        auto ptRec = candidate.pt();
        auto yRec = HfHelper::yD0(candidate);
        if (candidate.isRecoHfFlag() >= selectionFlagHf) {
          registry.fill(HIST("hPtVsYRecSigRecoHFFlag"), ptRec, yRec);
          registry.fill(HIST("hPtGenVsPtRecSig"), ptGen, ptRec);
          registry.fill(HIST("hYGenVsYRecSig"), yGen, yRec);
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hMassVsPtGenVsPtRecSig"), massD0, ptGen, ptRec);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hMassVsPtGenVsPtRecSig"), massD0bar, ptGen, ptRec);
          }
        }
        if (candidate.isRecoTopol() >= selectionTopol) {
          registry.fill(HIST("hPtVsYRecSigRecoTopol"), ptRec, yRec);
        }
        if (candidate.isRecoCand() >= selectionCand) {
          registry.fill(HIST("hPtVsYRecSigRecoCand"), ptRec, yRec);
        }
        if (candidate.isRecoPid() >= selectionPid) {
          registry.fill(HIST("hPtVsYRecSig_RecoPID"), ptRec, yRec);
        }
        if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0bar) {
          registry.fill(HIST("hPtVsYRecSigReco"), ptRec, yRec); // rec. level pT
          registry.fill(HIST("hPtRecSig"), ptRec);
        }
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          if (candidate.isRecoHfFlag() >= selectionFlagHf) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoHFFlag"), ptRec, yRec);
          }
          if (candidate.isRecoTopol() >= selectionTopol) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoTopol"), ptRec, yRec);
          }
          if (candidate.isRecoCand() >= selectionCand) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoCand"), ptRec, yRec);
          }
          if (candidate.isRecoPid() >= selectionPid) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoPID"), ptRec, yRec);
          }
          if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hPtVsYRecSigPromptReco"), ptRec, yRec); // rec. level pT, prompt
            registry.fill(HIST("hPtRecSigPrompt"), ptRec);
          }
        } else {
          if (candidate.isRecoHfFlag() >= selectionFlagHf) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoHFFlag"), ptRec, yRec);
          }
          if (candidate.isRecoTopol() >= selectionTopol) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoTopol"), ptRec, yRec);
          }
          if (candidate.isRecoCand() >= selectionCand) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoCand"), ptRec, yRec);
          }
          if (candidate.isRecoPid() >= selectionPid) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoPID"), ptRec, yRec);
          }
          if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hPtVsYRecSigNonPromptReco"), ptRec, yRec); // rec. level pT, non-prompt
            registry.fill(HIST("hPtRecSigNonPrompt"), ptRec);
          }
        }
        registry.fill(HIST("hCPARecSig"), candidate.cpa());
        registry.fill(HIST("hEtaRecSig"), candidate.eta());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa());
        registry.fill(HIST("hEtaRecBg"), candidate.eta());
      }
      auto ptCandidate = candidate.pt();
      auto ptProng0 = candidate.ptProng0();
      auto ptProng1 = candidate.ptProng1();
      auto rapidityCandidate = HfHelper::yD0(candidate);
      auto declengthCandidate = candidate.decayLength();
      auto declengthxyCandidate = candidate.decayLengthXY();
      auto normaliseddeclengthCandidate = candidate.decayLengthNormalised();
      auto normaliseddeclengthxyCandidate = candidate.decayLengthXYNormalised();
      auto d0Prong0 = candidate.impactParameter0();
      auto d0Prong1 = candidate.impactParameter1();
      auto d0d0Candidate = candidate.impactParameterProduct();
      auto ctsCandidate = HfHelper::cosThetaStarD0(candidate);
      auto ctCandidate = HfHelper::ctD0(candidate);
      auto cpaCandidate = candidate.cpa();
      auto cpaxyCandidate = candidate.cpaXY();
      int const minItsClustersOfProngs = std::min(trackPos.itsNCls(), trackNeg.itsNCls());
      int const minTpcCrossedRowsOfProngs = std::min(trackPos.tpcNClsCrossedRows(), trackNeg.tpcNClsCrossedRows());
      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMassSigBkgD0"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          registry.fill(HIST("hPtProng0Sig"), ptProng0, rapidityCandidate);
          registry.fill(HIST("hPtProng1Sig"), ptProng1, rapidityCandidate);
          registry.fill(HIST("hDecLengthSig"), declengthCandidate, rapidityCandidate);
          registry.fill(HIST("hDecLengthXYSig"), declengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthSig"), normaliseddeclengthCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthXYSig"), normaliseddeclengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hd0Prong0Sig"), d0Prong0, rapidityCandidate);
          registry.fill(HIST("hd0Prong1Sig"), d0Prong1, rapidityCandidate);
          registry.fill(HIST("hd0d0Sig"), d0d0Candidate, rapidityCandidate);
          registry.fill(HIST("hCTSSig"), ctsCandidate, rapidityCandidate);
          registry.fill(HIST("hCtSig"), ctCandidate, rapidityCandidate);
          registry.fill(HIST("hCPASig"), cpaCandidate, rapidityCandidate);
          registry.fill(HIST("hCPAxySig"), cpaxyCandidate, rapidityCandidate);
          registry.fill(HIST("hd0Prong0VsPtSig"), d0Prong0, ptCandidate);
          registry.fill(HIST("hd0Prong1VsPtSig"), d0Prong1, ptCandidate);
          registry.fill(HIST("hd0d0VsPtSig"), d0d0Candidate, ptCandidate);
          registry.fill(HIST("hCTSVsPtSig"), ctsCandidate, ptCandidate);
          registry.fill(HIST("hCPAVsPtSig"), cpaCandidate, ptCandidate);
          registry.fill(HIST("hCPAXYVsPtSig"), cpaxyCandidate, ptCandidate);
          registry.fill(HIST("hNormalisedDecLengthxyVsPtSig"), normaliseddeclengthxyCandidate, ptCandidate);
          registry.fill(HIST("hDecLengthVsPtSig"), declengthCandidate, ptCandidate);
          registry.fill(HIST("hDecLengthxyVsPtSig"), declengthxyCandidate, ptCandidate);
          registry.fill(HIST("hMassSigD0"), massD0, ptCandidate, rapidityCandidate);
          if constexpr (ApplyMl) {
            if (storeCentrality && storeOccupancyAndIR) {
              registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent, occ, ir);
            } else if (storeCentrality && !storeOccupancyAndIR) {
              registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent);
            } else if (!storeCentrality && storeOccupancyAndIR) {
              registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, occ, ir);
            } else if (storeTrackQuality) {
              registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
            } else {
              registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors);
            }
          } else {
            if (storeCentrality && storeOccupancyAndIR) {
              registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent, occ, ir);
            } else if (storeCentrality && !storeOccupancyAndIR) {
              registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent);
            } else if (!storeCentrality && storeOccupancyAndIR) {
              registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, occ, ir);
            } else if (storeTrackQuality) {
              registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
            } else {
              registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors);
            }
          }
        } else {
          registry.fill(HIST("hPtProng0Bkg"), ptProng0, rapidityCandidate);
          registry.fill(HIST("hPtProng1Bkg"), ptProng1, rapidityCandidate);
          registry.fill(HIST("hDecLengthBkg"), declengthCandidate, rapidityCandidate);
          registry.fill(HIST("hDecLengthXYBkg"), declengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthBkg"), normaliseddeclengthCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthXYBkg"), normaliseddeclengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hd0Prong0Bkg"), d0Prong0, rapidityCandidate);
          registry.fill(HIST("hd0Prong1Bkg"), d0Prong1, rapidityCandidate);
          registry.fill(HIST("hd0d0Bkg"), d0d0Candidate, rapidityCandidate);
          registry.fill(HIST("hCTSBkg"), ctsCandidate, rapidityCandidate);
          registry.fill(HIST("hCtBkg"), ctCandidate, rapidityCandidate);
          registry.fill(HIST("hCPABkg"), cpaCandidate, rapidityCandidate);
          registry.fill(HIST("hCPAxyBkg"), cpaxyCandidate, rapidityCandidate);
          registry.fill(HIST("hMassBkgD0"), massD0, ptCandidate, rapidityCandidate);
          if (candidate.flagMcMatchRec() == -o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
            registry.fill(HIST("hMassReflBkgD0"), massD0, ptCandidate, rapidityCandidate);
            if constexpr (ApplyMl) {
              if (storeCentrality && storeOccupancyAndIR) {
                registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent, occ, ir);
              } else if (storeCentrality && !storeOccupancyAndIR) {
                registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent);
              } else if (!storeCentrality && storeOccupancyAndIR) {
                registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, occ, ir);
              } else if (storeTrackQuality) {
                registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
              } else {
                registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0()[0], candidate.mlProbD0()[1], candidate.mlProbD0()[2], massD0, ptCandidate, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors);
              }
            } else {
              if (storeCentrality && storeOccupancyAndIR) {
                registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent, occ, ir);
              } else if (storeCentrality && !storeOccupancyAndIR) {
                registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent);
              } else if (!storeCentrality && storeOccupancyAndIR) {
                registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, occ, ir);
              } else if (storeTrackQuality) {
                registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
              } else {
                registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0, ptCandidate, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors);
              }
            }
          }
        }
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMassSigBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == -o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          registry.fill(HIST("hMassSigD0bar"), massD0bar, ptCandidate, rapidityCandidate);
          if constexpr (ApplyMl) {
            if (storeCentrality && storeOccupancyAndIR) {
              registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent, occ, ir);
            } else if (storeCentrality && !storeOccupancyAndIR) {
              registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent);
            } else if (!storeCentrality && storeOccupancyAndIR) {
              registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, occ, ir);
            } else if (storeTrackQuality) {
              registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
            } else {
              registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors);
            }
          } else {
            if (storeCentrality && storeOccupancyAndIR) {
              registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent, occ, ir);
            } else if (storeCentrality && !storeOccupancyAndIR) {
              registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent);
            } else if (!storeCentrality && storeOccupancyAndIR) {
              registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, occ, ir);
            } else if (storeTrackQuality) {
              registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
            } else {
              registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors);
            }
          }
        } else {
          registry.fill(HIST("hMassBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
          if (candidate.flagMcMatchRec() == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
            registry.fill(HIST("hMassReflBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
            if constexpr (ApplyMl) {
              if (storeCentrality && storeOccupancyAndIR) {
                registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent, occ, ir);
              } else if (storeCentrality && !storeOccupancyAndIR) {
                registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent);
              } else if (!storeCentrality && storeOccupancyAndIR) {
                registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, occ, ir);
              } else if (storeTrackQuality) {
                registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
              } else {
                registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"), candidate.mlProbD0bar()[0], candidate.mlProbD0bar()[1], candidate.mlProbD0bar()[2], massD0bar, ptCandidate, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors);
              }
            } else {
              if (storeCentrality && storeOccupancyAndIR) {
                registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent, occ, ir);
              } else if (storeCentrality && !storeOccupancyAndIR) {
                registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, cent);
              } else if (!storeCentrality && storeOccupancyAndIR) {
                registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, occ, ir);
              } else if (storeTrackQuality) {
                registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors, minItsClustersOfProngs, minTpcCrossedRowsOfProngs);
              } else {
                registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"), massD0bar, ptCandidate, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec(), numPvContributors);
              }
            }
          }
        }
      }
    }
    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
        if (yCandGenMax >= 0. && std::abs(RecoDecay::y(particle.pVector(), o2::constants::physics::MassD0)) > yCandGenMax) {
          continue;
        }
        float ptGenB = -1;
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassD0);
        registry.fill(HIST("hPtGen"), ptGen);
        registry.fill(HIST("hPtVsYGen"), ptGen, yGen);

        unsigned maxNumContrib = 0;
        float cent{-1.f};
        float occ{-1.f};
        if constexpr (std::is_same_v<CollType, CollisionsWithMcLabelsCent>) {
          const auto& recoCollsPerMcCollCent = collisions.sliceBy(colPerMcCollisionCent, particle.mcCollision().globalIndex());
          for (const auto& recCol : recoCollsPerMcCollCent) {
            maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
          }
          if (storeCentrality && centEstimator != CentralityEstimator::None) {
            cent = getCentralityGenColl(recoCollsPerMcCollCent, centEstimator);
          }
          if (storeOccupancyAndIR && occEstimator != OccupancyEstimator::None) {
            occ = o2::hf_occupancy::getOccupancyGenColl(recoCollsPerMcCollCent, occEstimator);
          }
        } else {
          const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, particle.mcCollision().globalIndex());
          for (const auto& recCol : recoCollsPerMcColl) {
            maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
          }
        }

        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hPtGenPrompt"), ptGen);
          registry.fill(HIST("hYGenPrompt"), yGen);
          registry.fill(HIST("hPtVsYGenPrompt"), ptGen, yGen);
          if (storeCentrality && storeOccupancyAndIR) {
            registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, 1, maxNumContrib, cent, occ);
          } else if (storeCentrality && !storeOccupancyAndIR) {
            registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, 1, maxNumContrib, cent);
          } else if (!storeCentrality && storeOccupancyAndIR) {
            registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, 1, maxNumContrib, occ);
          } else {
            registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, 1, maxNumContrib);
          }
        } else {
          ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
          registry.fill(HIST("hPtGenNonPrompt"), ptGen);
          registry.fill(HIST("hYGenNonPrompt"), yGen);
          registry.fill(HIST("hPtVsYGenNonPrompt"), ptGen, yGen);
          if (storeCentrality && storeOccupancyAndIR) {
            registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, 2, maxNumContrib, cent, occ);
          } else if (storeCentrality && !storeOccupancyAndIR) {
            registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, 2, maxNumContrib, cent);
          } else if (!storeCentrality && storeOccupancyAndIR) {
            registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, 2, maxNumContrib, occ);
          } else {
            registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, 2, maxNumContrib);
          }
        }
        registry.fill(HIST("hEtaGen"), particle.eta());
      }
    }
  }

  void processMcWithDCAFitterN(D0CandidatesMc const&,
                               soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
                               TracksSelQuality const& tracks,
                               CollisionsWithMcLabels const& collisions,
                               aod::McCollisions const& mcCollisions,
                               aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, false>(selectedD0CandidatesMc, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processMcWithDCAFitterN, "Process MC with DCAFitterN", false);

  void processMcWithDCAFitterNCent(D0CandidatesMc const&,
                                   soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
                                   TracksSelQuality const& tracks,
                                   CollisionsWithMcLabelsCent const& collisions,
                                   aod::McCollisions const& mcCollisions,
                                   aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, false>(selectedD0CandidatesMc, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processMcWithDCAFitterNCent, "Process MC with DCAFitterN and centrality", false);

  void processMcWithKFParticle(D0CandidatesMcKF const&,
                               soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
                               TracksSelQuality const& tracks,
                               CollisionsWithMcLabels const& collisions,
                               aod::McCollisions const& mcCollisions,
                               aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, false>(selectedD0CandidatesMcKF, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processMcWithKFParticle, "Process MC with KFParticle", false);
  // TODO: add the processMcWithKFParticleCent

  void processMcWithDCAFitterNMl(D0CandidatesMlMc const&,
                                 soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
                                 TracksSelQuality const& tracks,
                                 CollisionsWithMcLabels const& collisions,
                                 aod::McCollisions const& mcCollisions,
                                 aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, true>(selectedD0CandidatesMlMc, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processMcWithDCAFitterNMl, "Process MC with DCAFitterN and ML selection", false);

  void processMcWithDCAFitterNMlCent(D0CandidatesMlMc const&,
                                     soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
                                     TracksSelQuality const& tracks,
                                     CollisionsWithMcLabelsCent const& collisions,
                                     aod::McCollisions const& mcCollisions,
                                     aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, true>(selectedD0CandidatesMlMc, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processMcWithDCAFitterNMlCent, "Process MC with DCAFitterN and ML selection and centrality", false);

  void processMcWithKFParticleMl(D0CandidatesMlMcKF const&,
                                 soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
                                 TracksSelQuality const& tracks,
                                 CollisionsWithMcLabels const& collisions,
                                 aod::McCollisions const& mcCollisions,
                                 aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, true>(selectedD0CandidatesMlMcKF, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0, processMcWithKFParticleMl, "Process MC with KFParticle and ML selections", false);
  // TODO: add the processMcWithKFParticleMlCent
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskD0>(cfgc)};
}
