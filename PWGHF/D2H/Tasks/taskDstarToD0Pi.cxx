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

/// \file taskDstarToD0Pi.cxx
/// \file D* analysis task

/// \author Deependra Sharma <deependra.sharma@cern.ch>, IITB
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

/// \brief Dstar production analysis task (With and Without ML)

#include <algorithm>
#include <utility>
#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

/// Dstar analysis task
struct HfTaskDstarToD0Pi {
  Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 & Pi"};
  Configurable<double> yCandDstarRecoMax{"yCandDstarRecoMax", 0.8, "max. candidate Dstar rapidity"};
  Configurable<double> yCandDstarGenMax{"yCandDstarGenMax", 0.5, "max. rapidity of Generator level Particle"};
  Configurable<bool> selectionFlagHfD0ToPiK{"selectionFlagHfD0ToPiK", true, "Selection Flag for HF flagged candidates"};
  Configurable<std::vector<double>> ptBins{"ptBins", std::vector<double>{hf_cuts_dstar_to_d0_pi::vecBinsPt}, "pT bin limits for Dstar"};

  SliceCache cache;

  using CandDstarWSelFlag = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi>;
  using CandDstarWSelFlagWMl = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi>;
  /// @brief specially for MC data
  // full reconstructed Dstar candidate
  using CandDstarWSelFlagMcRec = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfCandDstarMcRec>;
  using CandDstarWSelFlagWMlMcRec = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfCandDstarMcRec, aod::HfMlDstarToD0Pi>;
  using CandDstarMcGen = soa::Join<aod::McParticles, aod::HfCandDstarMcGen>;

  using CollisionsWCent = soa::Join<aod::Collisions, aod::CentFT0Ms>;
  using CollisionsWCentMcLabel = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Ms>;

  Filter candFilter = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;
  Preslice<soa::Filtered<CandDstarWSelFlag>> preslicSelectedCandDstarPerCol = aod::hf_cand::collisionId;
  Preslice<soa::Filtered<CandDstarWSelFlagWMl>> preslicSelectedCandDstarPerColWMl = aod::hf_cand::collisionId;

  PresliceUnsorted<CollisionsWCentMcLabel> colsPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  Partition<CandDstarWSelFlagMcRec> rowsSelectedCandDstarMcRec = aod::hf_sel_candidate_dstar::isRecoD0Flag == selectionFlagHfD0ToPiK;
  Partition<CandDstarWSelFlagWMlMcRec> rowsSelectedCandDstarMcRecWMl = aod::hf_sel_candidate_dstar::isRecoD0Flag == selectionFlagHfD0ToPiK;

  ConfigurableAxis binningImpactParam{"binningImpactParam", {1000, 0.1, -0.1}, " Bins of Impact Parameter"};
  ConfigurableAxis binningDecayLength{"binningDecayLength", {1000, 0.0, 0.7}, "Bins of Decay Length"};
  ConfigurableAxis binningNormDecayLength{"binningNormDecayLength", {1000, 0.0, 40.0}, "Bins of Normalised Decay Length"};
  ConfigurableAxis binningCentrality{"binningCentrality", {VARIABLE_WIDTH, 0.0, 1.0, 10.0, 30.0, 50.0, 70.0, 100.0}, "centrality binning"};
  ConfigurableAxis binningDeltaInvMass{"binningDeltaInvMass", {100, 0.13, 0.16}, "Bins of Delta InvMass of Dstar"};
  ConfigurableAxis binningBkgBDTScore{"binningBkgBDTScore", {100, 0.0f, 1.0f}, "Bins for background BDT Score"};
  ConfigurableAxis binningSigBDTScore{"binningSigBDTScore", {100, 0.0f, 1.0f}, "Bins for Signal (Prompts + Non Prompt) BDT Score"};

  HistogramRegistry registry{
    "registry",
    {{"QA/hPtDstar", "Dstar Candidates; Dstar candidate #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtD0", "D0 Candiades; D0 Candidate #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtProng0D0", "Prong0 of D0 Candidates; Prong0 #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtProng1D0", "Prong1 of D0 Candidates; Prong1 #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtProng0D0Bar", "Prong0 of D0Bar Candidates; Prong0 #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtProng1D0Bar", "Prong1 of D0Bar Candidates; Prong1 #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"QA/hPtSoftPi", "Soft Pi ; Soft Pi #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, 0., 1.}}}},
     {"QA/hDeltaCentGen", "#{Delta}Cent % distribution of Collisions having same MC Collision;FT0M #{Delta}Cent %; Counts", {HistType::kTH1F, {{100, 0.0, 100.0}}}}}};

  void init(InitContext&)
  {
    if ((doprocessDataWoML && doprocessDataWML) || (doprocessMcWoMl && doprocessMcWML) || (doprocessDataWoML && doprocessMcWML) || (doprocessDataWML && doprocessMcWoMl)) {
      LOGP(fatal, "Only Without-ML or With-ML functions should be enabled at a time! Please check your configuration!");
    }
    auto vecPtBins = (std::vector<double>)ptBins;

    AxisSpec axisImpactParam = {binningImpactParam, "impact parameter (cm)"};
    AxisSpec axisDecayLength = {binningDecayLength, " decay length (cm)"};
    AxisSpec axisNormDecayLength = {binningNormDecayLength, "normalised decay length (cm)"};
    AxisSpec axisCentrality = {binningCentrality, "centrality (%)"};
    AxisSpec axisDeltaInvMass = {binningDeltaInvMass, "#Delta #it{M}_{inv} D*"};
    AxisSpec axisBDTScorePrompt = {binningSigBDTScore, "BDT Score for Prompt Cand"};
    AxisSpec axisBDTScoreNonPrompt = {binningSigBDTScore, "BDT Score for Non-Prompt Cand"};
    AxisSpec axisBDTScoreBackground = {binningBkgBDTScore, "BDT Score for Background Cand"};

    registry.add("Yield/hDeltaInvMassDstar3D", "#Delta #it{M}_{inv} D* Candidate; inv. mass ((#pi #pi k) - (#pi k)) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c}); FT0M centrality", {HistType::kTH3F, {{axisDeltaInvMass}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {axisCentrality}}}, true);
    registry.add("Yield/hDeltaInvMassDstar2D", "#Delta #it{M}_{inv} D* Candidate; inv. mass ((#pi #pi k) - (#pi k)) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{axisDeltaInvMass}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}}, true);
    registry.add("Yield/hDeltaInvMassDstar1D", "#Delta #it{M}_{inv} D* Candidate; inv. mass ((#pi #pi k) - (#pi k)) (GeV/#it{c}^{2}); entries", {HistType::kTH1F, {{axisDeltaInvMass}}}, true);
    registry.add("Yield/hInvMassDstar", "#Delta #it{M}_{inv} D* Candidate; inv. mass (#pi #pi k) (GeV/#it{c}^{2}); entries", {HistType::kTH1F, {{500, 0., 5.0}}}, true);
    registry.add("Yield/hInvMassD0", "#it{M}_{inv}D^{0} candidate;#it{M}_{inv} D^{0} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{500, 0., 5.0}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}}, true);
    // only QA
    registry.add("QA/hEtaDstar", "D* Candidate; D* Candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hCtD0", "D0 Candidate; D0 Candidate's proper life time #it{c}#tau (cm) ; entries ", {HistType::kTH2F, {{1000, -0.1, 14.}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthD0", "D0 Candidate; D0 Candidate's decay length (cm); entries", {HistType::kTH2F, {axisDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthXYD0", "D0 Candidate; D0 Candidate's decay length xy (cm);  entries", {HistType::kTH2F, {axisDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthNormalisedD0", "D0 Candidates;Do Candidate's normalised decay length (cm); entries", {HistType::kTH2F, {axisNormDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthXYNormalisedD0", "D0 candidate; D0 Candidate's normalised decay length xy (cm); entries", {HistType::kTH2F, {axisNormDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hCPAD0", "D0 Candidates; D0 Candidate's cosine pointing angle; entries", {HistType::kTH2F, {{110, -1., 1.}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hCPAxyD0", "D0 candidates; D0 Candidate's cosine of pointing angle xy; entries", {HistType::kTH2F, {{110, -1., 1.}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hImpactParameterXYD0", "D0 Candidates; D0 Candidate's reconstructed impact parameter xy (cm); entries", {HistType::kTH2F, {axisImpactParam, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDeltaIPMaxNormalisedD0", "D0 Candidate; D0 Candidate's Norm. Delta IP; entries", {HistType::kTH2F, {{1000, -20., 20.}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hSqSumProngsImpactParameterD0", "D0 Candidates; Sqr Sum of Impact params of D0 Prongs; entries ", {HistType::kTH2F, {{1000, 0., 0.25}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthErrorD0", "D0 Candidates; D0 Candidate's Decay Length Error (cm); entries", {HistType::kTH2F, {axisDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hDecayLengthXYErrorD0", "D0 Candidates; D0 Candidate's Decay Length Error XY (cm); entries", {HistType::kTH2F, {axisDecayLength, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hImpactParameterError", "D0 Prongs; Impact param error of different D0 Prongs (cm); entries", {HistType::kTH2F, {axisImpactParam, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hd0Prong0", "Prong0; DCAxy of Prong0 (cm); entries", {HistType::kTH2F, {axisImpactParam, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hd0Prong1", "Prong1; DCAxy of Prong1 (cm); entries", {HistType::kTH2F, {axisImpactParam, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hd0ProngSoftPi", "ProngSoftPi; DCAxy of Prong Soft Pi (cm); entries", {HistType::kTH2F, {axisImpactParam, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    // MC Matching at Reconstruction Level Successful
    registry.add("QA/hCPASkimD0RecSig", "MC Matched Skimed D* Candidates at Reconstruction Level; cosine of pointing angle", {HistType::kTH1F, {{110, -1.1, 1.1}}});
    registry.add("QA/hEtaSkimD0RecSig", "MC Matched Skimed D* Candidates at Reconstruction Level; #it{#eta} of D0 Prong", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("QA/hEtaSkimDstarRecSig", "MC Matched Skimed D* Candidates at Reconstruction Level; #it{#eta} of D* Candidate", {HistType::kTH1F, {{100, -2., 2.}}});
    // pt vs y
    registry.add("QA/hPtSkimDstarGenSig", "MC Matched Skimed D* Reconstructed Candidates at Generator Level; #it{p}_{T} of D* at Generator Level (GeV/#it{c})", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}}, true);
    registry.add("Efficiency/hPtVsCentSkimDstarGenSig", "MC Matched Skimed D* Reconstructed Candidates at Generator Level; #it{p}_{T} of D* at Generator Level (GeV/#it{c}); Centrality (%)", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {axisCentrality}}}, true);
    registry.add("QA/hPtVsYSkimDstarRecSig", "MC Matched Skimed D* Candidates at Reconstruction Level; #it{p}_{T} of D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoTopolDstarRecSig", "MC Matched RecoTopol D* Candidates at Reconstruction Level; #it{p}_{T} of D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoPidDstarRecSig", "MC Matched RecoPid D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtFullRecoDstarRecSig", "MC Matched FullReco D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c})", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Efficiency/hPtVsCentFullRecoDstarRecSig", "MC Matched  FullReco D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); Centrality (%)", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {axisCentrality}}}, true);
    // Only Prompt RecSig
    registry.add("QA/hPtVsYSkimPromptDstarRecSig", "MC Matched Skimed Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}; #it{y})", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoTopolPromptDstarRecSig", "MC Matched RecoTopol Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoPidPromptDstarRecSig", "MC Matched RecoPid Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtFullRecoPromptDstarRecSig", "MC Matched FullReco Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c})", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    // Only Non-Prompt RecSig
    registry.add("QA/hPtVsYSkimNonPromptDstarRecSig", "MC Matched Skimed Non-Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoTopolNonPromptDstarRecSig", "MC Matched RecoTopol Non-Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtVsYRecoPidNonPromptDstarRecSig", "MC Matched RecoPid Non-Prompt D*  Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c}); #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("QA/hPtFullRecoNonPromptDstarRecSig", "MC Matched FullReco Non-Prompt D* Candidates at Reconstruction Level; #it{p}_{T} of  D* at Reconstruction Level (GeV/#it{c})", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    // MC Matching UnSuccessful
    registry.add("QA/hCPASkimD0RecBg", "MC UnMatched Skimmed D0 Candidates at Reconstruction Level; cosine of pointing angle", {HistType::kTH1F, {{110, -1., 1.}}});
    registry.add("QA/hEtaSkimD0RecBg", "MC UnMatched Skimmed D0 Candidates at Reconstruction Level; #it{#eta} of D0 Prong", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("QA/hEtaSkimDstarRecBg", "MC UnMatched Skimmed D* Candidates at Reconstruction Level; #it{#eta} of D* Candidate", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("QA/hPtSkimD0RecBg", "MC UnMatched Skimmed D0 Candidates at Reconstruction Level;  #it{p}_{T} of  D0", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hPtSkimDstarRecBg", "MC UnMatched Skimmed D* Candidates at Reconstruction Level; #it{p}_{T} of  D*", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    // MC Matching at Generator level Successful
    registry.add("QA/hEtaDstarGen", "MC Matched D* Candidates at Generator Level; #it{#eta}", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("QA/hPtDstarGen", "MC Matched D* Candidates at Generator Level; #it{p}_{T} of  D*", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Efficiency/hPtVsCentDstarGen", "MC Matched D* Candidates at Generator Level; #it{p}_{T} of  D*;Centrality (%) ", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {axisCentrality}}}, true);
    registry.add("QA/hPtVsYDstarGen", "MC Matched D* Candidates at Generator Level; #it{p}_{T} of  D*; #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    // Prompt Gen
    registry.add("QA/hPtPromptDstarGen", "MC Matched Prompt D* Candidates at Generator Level; #it{p}_{T} of  D*", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hPtVsYPromptDstarGen", "MC Matched Prompt D* Candidates at Generator Level; #it{p}_{T} of  D*; #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    // Non Prmpt Gen
    registry.add("QA/hPtNonPromptDstarGen", "MC Matched Non-Prompt D* Candidates at Generator Level; #it{p}_{T} of  D*", {HistType::kTH1F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("QA/hPtVsYNonPromptDstarGen", "MC Matched Non-Prompt D* Candidates at Generator Level; #it{p}_{T} of  D*; #it{y}", {HistType::kTH2F, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});

    // Checking PV contributors from Data as well MC rec
    registry.add("Efficiency/hNumPvContributorsAll", "PV Contributors; PV Contributor; FT0M Centrality", {HistType::kTH2F, {{100, 0, 300}, {axisCentrality}}}, true);
    registry.add("Efficiency/hNumPvContributorsCand", "PV Contributors; PV Contributor; FT0M Centrality", {HistType::kTH2F, {{100, 0, 300}, {axisCentrality}}}, true);
    registry.add("Efficiency/hNumPvContributorsCandInMass", "PV Contributors; PV Contributor; FT0M Centrality", {HistType::kTH2F, {{100, 0, 300}, {axisCentrality}}}, true);

    // BDT Score (axisBDTScoreBackground, axisBDTScorePrompt, axisBDTScoreNonPrompt)
    if (doprocessDataWML) {
      registry.add("Yield/hDeltaInvMassVsPtVsCentVsBDTScore", "#Delta #it{M}_{inv} Vs Pt Vs Cent Vs BDTScore", {HistType::kTHnSparseF, {{axisDeltaInvMass}, {vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {axisCentrality}, {axisBDTScoreBackground}, {axisBDTScorePrompt}, {axisBDTScoreNonPrompt}}});
    }
    if (doprocessMcWML) {
      registry.add("Efficiency/hPtVsCentVsBDTScore", "Pt Vs Cent Vs BDTScore", {HistType::kTHnSparseF, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {axisCentrality}, {axisBDTScoreBackground}, {axisBDTScorePrompt}, {axisBDTScoreNonPrompt}}});
      registry.add("Efficiency/hPtPromptVsCentVsBDTScore", "Pt Vs Cent Vs BDTScore", {HistType::kTHnSparseF, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {axisCentrality}, {axisBDTScoreBackground}, {axisBDTScorePrompt}, {axisBDTScoreNonPrompt}}});
      registry.add("Efficiency/hPtNonPromptVsCentVsBDTScore", "Pt Vs Cent Vs BDTScore", {HistType::kTHnSparseF, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {axisCentrality}, {axisBDTScoreBackground}, {axisBDTScorePrompt}, {axisBDTScoreNonPrompt}}});
      // registry.add("Efficiency/hPtBkgVsCentVsBDTScore", "Pt Vs Cent Vs BDTScore", {HistType::kTHnSparseF, {{vecPtBins, "#it{p}_{T} (GeV/#it{c})"}, {axisCentrality}, {axisBDTScoreBackground}, {axisBDTScorePrompt}, {axisBDTScoreNonPrompt}}});
    }
  }

  // Comparator function to sort based on the second argument of a tuple
  static bool compare(const std::pair<soa::Filtered<CollisionsWCentMcLabel>::iterator, int>& a, const std::pair<soa::Filtered<CollisionsWCentMcLabel>::iterator, int>& b)
  {
    return a.second > b.second;
  }

  /// @brief This function runs over Data to obatin yield
  /// @tparam T1 type of the candidate
  /// @tparam T2 type of preslice used to slice the candidate table
  /// @tparam applyMl  a boolean to apply ML or not
  /// @param cols reconstructed collision with centrality
  /// @param selectedCands selected candidates with selection flag
  /// @param preslice preslice to slice
  template <bool applyMl, typename T1, typename T2>
  void runTaskDstar(CollisionsWCent const& cols, T1 selectedCands, T2 preslice)
  {
    for (const auto& col : cols) {
      auto nPVContributors = col.numContrib();
      auto centrality = col.centFT0M();
      registry.fill(HIST("Efficiency/hNumPvContributorsAll"), nPVContributors, centrality);

      auto gIndexCol = col.globalIndex();
      auto selectedCandsCurrentCol = selectedCands.sliceBy(preslice, gIndexCol);
      auto nCandsCurrentCol = selectedCandsCurrentCol.size();

      if (nCandsCurrentCol > 0) {
        // LOGF(debug, "size of selectedCandsCurrentCol: %d", nCandsCurrentCol);
        registry.fill(HIST("Efficiency/hNumPvContributorsCand"), nPVContributors, centrality);
      }

      int nCandsSignalRegion = 0;
      for (const auto& candDstar : selectedCandsCurrentCol) {
        auto yDstar = candDstar.y(constants::physics::MassDStar);
        if (yCandDstarRecoMax >= 0. && std::abs(yDstar) > yCandDstarRecoMax) {
          continue;
        }

        registry.fill(HIST("QA/hPtDstar"), candDstar.pt());
        registry.fill(HIST("QA/hPtD0"), candDstar.ptD0());
        registry.fill(HIST("QA/hPtSoftPi"), candDstar.ptSoftPi());
        registry.fill(HIST("QA/hEtaDstar"), candDstar.eta(), candDstar.pt());
        registry.fill(HIST("QA/hCtD0"), candDstar.ctD0(), candDstar.pt());
        registry.fill(HIST("QA/hDecayLengthD0"), candDstar.decayLengthD0(), candDstar.pt());
        registry.fill(HIST("QA/hDecayLengthXYD0"), candDstar.decayLengthXYD0(), candDstar.pt());
        registry.fill(HIST("QA/hDecayLengthNormalisedD0"), candDstar.decayLengthNormalisedD0(), candDstar.pt());
        registry.fill(HIST("QA/hDecayLengthXYNormalisedD0"), candDstar.decayLengthXYNormalisedD0(), candDstar.pt());
        registry.fill(HIST("QA/hCPAD0"), candDstar.cpaD0(), candDstar.pt());
        registry.fill(HIST("QA/hCPAxyD0"), candDstar.cpaXYD0(), candDstar.pt());
        registry.fill(HIST("QA/hImpactParameterXYD0"), candDstar.impactParameterXYD0(), candDstar.pt());
        registry.fill(HIST("QA/hDeltaIPMaxNormalisedD0"), candDstar.deltaIPNormalisedMaxD0(), candDstar.pt());
        registry.fill(HIST("QA/hSqSumProngsImpactParameterD0"), candDstar.impactParameterProngSqSumD0(), candDstar.pt());
        registry.fill(HIST("QA/hDecayLengthErrorD0"), candDstar.errorDecayLengthD0(), candDstar.pt());
        registry.fill(HIST("QA/hDecayLengthXYErrorD0"), candDstar.errorDecayLengthXYD0(), candDstar.pt());
        registry.fill(HIST("QA/hImpactParameterError"), candDstar.errorImpactParameter0(), candDstar.pt());
        registry.fill(HIST("QA/hImpactParameterError"), candDstar.errorImpactParameter1(), candDstar.pt());
        registry.fill(HIST("QA/hImpactParameterError"), candDstar.errorImpParamSoftPi(), candDstar.pt());
        registry.fill(HIST("QA/hd0Prong0"), candDstar.impactParameter0(), candDstar.pt());
        registry.fill(HIST("QA/hd0Prong1"), candDstar.impactParameter1(), candDstar.pt());
        registry.fill(HIST("QA/hd0ProngSoftPi"), candDstar.impParamSoftPi(), candDstar.pt());

        auto invDstar = candDstar.invMassDstar();
        auto invAntiDstar = candDstar.invMassAntiDstar();
        auto invD0 = candDstar.invMassD0();
        auto invD0Bar = candDstar.invMassD0Bar();

        auto signDstar = candDstar.signSoftPi();
        if (signDstar > 0) {
          auto deltaMDstar = std::abs(invDstar - invD0);
          if (0.142f < deltaMDstar && deltaMDstar < 0.15f) {
            nCandsSignalRegion++;
          }

          if constexpr (applyMl) {
            auto mlBdtScore = candDstar.mlProbDstarToD0Pi();
            registry.fill(HIST("Yield/hDeltaInvMassVsPtVsCentVsBDTScore"), deltaMDstar, candDstar.pt(), centrality, mlBdtScore[0], mlBdtScore[1], mlBdtScore[2]);
          }

          registry.fill(HIST("Yield/hDeltaInvMassDstar3D"), deltaMDstar, candDstar.pt(), centrality);
          registry.fill(HIST("Yield/hDeltaInvMassDstar2D"), deltaMDstar, candDstar.pt());
          registry.fill(HIST("Yield/hInvMassD0"), invD0, candDstar.ptD0());
          registry.fill(HIST("Yield/hDeltaInvMassDstar1D"), deltaMDstar);
          registry.fill(HIST("Yield/hInvMassDstar"), invDstar);
          // filling pt of two pronges of D0
          registry.fill(HIST("QA/hPtProng0D0"), candDstar.ptProng0());
          registry.fill(HIST("QA/hPtProng1D0"), candDstar.ptProng1());
        } else if (signDstar < 0) {
          auto deltaMAntiDstar = std::abs(invAntiDstar - invD0Bar);
          if (0.142f < deltaMAntiDstar && deltaMAntiDstar < 0.15f) {
            nCandsSignalRegion++;
          }

          if constexpr (applyMl) {
            auto mlBdtScore = candDstar.mlProbDstarToD0Pi();
            registry.fill(HIST("Yield/hDeltaInvMassVsPtVsCentVsBDTScore"), deltaMAntiDstar, candDstar.pt(), centrality, mlBdtScore[0], mlBdtScore[1], mlBdtScore[2]);
          }

          registry.fill(HIST("Yield/hDeltaInvMassDstar3D"), deltaMAntiDstar, candDstar.pt(), centrality);
          registry.fill(HIST("Yield/hDeltaInvMassDstar2D"), deltaMAntiDstar, candDstar.pt());
          registry.fill(HIST("Yield/hInvMassD0"), invD0Bar, candDstar.ptD0());
          registry.fill(HIST("Yield/hDeltaInvMassDstar1D"), deltaMAntiDstar);
          registry.fill(HIST("Yield/hInvMassDstar"), invAntiDstar);
          // filling pt of two pronges of D0Bar
          registry.fill(HIST("QA/hPtProng0D0Bar"), candDstar.ptProng0());
          registry.fill(HIST("QA/hPtProng1D0Bar"), candDstar.ptProng1());
        }
      } // candidate loop for current collision ends

      if (nCandsSignalRegion > 0) {
        registry.fill(HIST("Efficiency/hNumPvContributorsCandInMass"), nPVContributors, centrality);
      }
    } // collision loop ends
  }

  /// @brief This function runs over MC at reco level to obatin efficiency
  /// @tparam T1 type of the candidate table
  /// @tparam applyMl a boolean to apply ML or not
  /// @param candsMcRecSel reconstructed candidates with selection flag
  /// @param rowsMcPartilces generated particles  table
  template <bool applyMl, typename T1>
  void runMcRecTaskDstar(T1 const& candsMcRecSel, CandDstarMcGen const& rowsMcPartilces)
  {
    // LOGF(info, "Running MC Rec Task Dstar");
    int8_t signDstar = 0;
    // MC at Reconstruction level
    for (const auto& candDstarMcRec : candsMcRecSel) {
      // LOGF(info, "MC Rec Dstar loop");
      auto ptDstarRecSig = candDstarMcRec.pt();
      auto yDstarRecSig = candDstarMcRec.y(constants::physics::MassDStar);
      if (yCandDstarRecoMax >= 0. && std::abs(yDstarRecSig) > yCandDstarRecoMax) {
        continue;
      }
      auto collision = candDstarMcRec.template collision_as<CollisionsWCentMcLabel>();
      auto centrality = collision.centFT0M();                                                               // 0-100%
      if (TESTBIT(std::abs(candDstarMcRec.flagMcMatchRec()), aod::hf_cand_dstar::DecayType::DstarToD0Pi)) { // if MC matching is successful at Reconstruction Level
        // LOGF(info, "MC Rec Dstar loop MC Matched");
        // get MC Mother particle
        auto prong0 = candDstarMcRec.template prong0_as<aod::TracksWMc>();
        auto indexMother = RecoDecay::getMother(rowsMcPartilces, prong0.template mcParticle_as<CandDstarMcGen>(), o2::constants::physics::Pdg::kDStar, true, &signDstar, 2);
        auto particleMother = rowsMcPartilces.rawIteratorAt(indexMother);  // What is difference between rawIterator() or iteratorAt() methods?
        registry.fill(HIST("QA/hPtSkimDstarGenSig"), particleMother.pt()); // generator level pt
        registry.fill(HIST("Efficiency/hPtVsCentSkimDstarGenSig"), particleMother.pt(), centrality);

        registry.fill(HIST("QA/hPtVsYSkimDstarRecSig"), ptDstarRecSig, yDstarRecSig); // Skimed at level of trackIndexSkimCreator
        if (candDstarMcRec.isRecoTopol()) {                                           // if Topological selection are passed
          registry.fill(HIST("QA/hPtVsYRecoTopolDstarRecSig"), ptDstarRecSig, yDstarRecSig);
        }
        if (candDstarMcRec.isRecoPid()) { // if PID selection is passed
          registry.fill(HIST("QA/hPtVsYRecoPidDstarRecSig"), ptDstarRecSig, yDstarRecSig);
        }
        if (candDstarMcRec.isSelDstarToD0Pi()) { // if all selection passed
          registry.fill(HIST("QA/hPtFullRecoDstarRecSig"), ptDstarRecSig);
          registry.fill(HIST("Efficiency/hPtVsCentFullRecoDstarRecSig"), ptDstarRecSig, centrality);
          if constexpr (applyMl) {
            auto bdtScore = candDstarMcRec.mlProbDstarToD0Pi();
            registry.fill(HIST("Efficiency/hPtVsCentVsBDTScore"), ptDstarRecSig, centrality, bdtScore[0], bdtScore[1], bdtScore[2]);
          }
        }
        registry.fill(HIST("QA/hCPASkimD0RecSig"), candDstarMcRec.cpaD0());
        registry.fill(HIST("QA/hEtaSkimD0RecSig"), candDstarMcRec.etaD0());
        registry.fill(HIST("QA/hEtaSkimDstarRecSig"), candDstarMcRec.eta());

        // only prompt signal at reconstruction level
        if (candDstarMcRec.originMcRec() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("QA/hPtVsYSkimPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig); // Skimed at level of trackIndexSkimCreator
          if (candDstarMcRec.isRecoTopol()) {                                                 // if Topological selection are passed
            registry.fill(HIST("QA/hPtVsYRecoTopolPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig);
          }
          if (candDstarMcRec.isRecoPid()) { // if PID selection is passed
            registry.fill(HIST("QA/hPtVsYRecoPidPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig);
          }
          if (candDstarMcRec.isSelDstarToD0Pi()) { // if all selection passed
            registry.fill(HIST("QA/hPtFullRecoPromptDstarRecSig"), ptDstarRecSig);
            if constexpr (applyMl) {
              auto bdtScore = candDstarMcRec.mlProbDstarToD0Pi();
              registry.fill(HIST("Efficiency/hPtPromptVsCentVsBDTScore"), ptDstarRecSig, centrality, bdtScore[0], bdtScore[1], bdtScore[2]);
            }
          }
        } else if (candDstarMcRec.originMcRec() == RecoDecay::OriginType::NonPrompt) { // only non-prompt signal at reconstruction level
          registry.fill(HIST("QA/hPtVsYSkimNonPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig);
          if (candDstarMcRec.isRecoTopol()) { // if Topological selection are passed
            registry.fill(HIST("QA/hPtVsYRecoTopolNonPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig);
          }
          if (candDstarMcRec.isRecoPid()) { // if PID selection is passed
            registry.fill(HIST("QA/hPtVsYRecoPidNonPromptDstarRecSig"), ptDstarRecSig, yDstarRecSig);
          }
          if (candDstarMcRec.isSelDstarToD0Pi()) { // if all selection passed
            registry.fill(HIST("QA/hPtFullRecoNonPromptDstarRecSig"), ptDstarRecSig);
            if constexpr (applyMl) {
              auto bdtScore = candDstarMcRec.mlProbDstarToD0Pi();
              registry.fill(HIST("Efficiency/hPtNonPromptVsCentVsBDTScore"), ptDstarRecSig, centrality, bdtScore[0], bdtScore[1], bdtScore[2]);
            }
          }
        }
      } else { // MC Unmatched (Baground at Reconstruction Level)
        registry.fill(HIST("QA/hCPASkimD0RecBg"), candDstarMcRec.cpaD0());
        registry.fill(HIST("QA/hEtaSkimD0RecBg"), candDstarMcRec.etaD0());
        registry.fill(HIST("QA/hEtaSkimDstarRecBg"), candDstarMcRec.eta());
        if (candDstarMcRec.isSelDstarToD0Pi()) {
          registry.fill(HIST("QA/hPtSkimDstarRecBg"), candDstarMcRec.pt());
        }
      }
    } // candidate loop ends
  }

  /// @brief This function runs over MC at gen level to obatin efficiency
  /// @param collisions reconstructed collision with centrality
  /// @param rowsMcPartilces generated particles  table
  void runMcGenTaskDstar(CollisionsWCentMcLabel const& collisions, CandDstarMcGen const& rowsMcPartilces)
  {
    // MC Gen level
    for (auto const& mcParticle : rowsMcPartilces) {
      if (TESTBIT(std::abs(mcParticle.flagMcMatchGen()), aod::hf_cand_dstar::DecayType::DstarToD0Pi)) { // MC Matching is successful at Generator Level
        auto ptGen = mcParticle.pt();
        auto yGen = RecoDecay::y(mcParticle.pVector(), o2::constants::physics::MassDStar);
        if (yCandDstarGenMax >= 0. && std::abs(yGen) > yCandDstarGenMax) {
          continue;
        }
        auto mcCollision = mcParticle.template mcCollision_as<aod::McCollisions>();
        auto recCollisions = collisions.sliceBy(colsPerMcCollision, mcCollision.globalIndex());
        // looking if a generated collision reconstructed more than a times.
        if (recCollisions.size() > 1) {
          for (const auto& [c1, c2] : combinations(CombinationsStrictlyUpperIndexPolicy(recCollisions, recCollisions))) {
            auto deltaCent = std::abs(c1.centFT0M() - c2.centFT0M());
            registry.fill(HIST("QA/hDeltaCentGen"), deltaCent);
          }
        }
        float centFT0MGen;
        // assigning centrality to MC Collision using max FT0M amplitute from Reconstructed collisions
        if (recCollisions.size()) {
          std::vector<std::pair<soa::Filtered<CollisionsWCentMcLabel>::iterator, int>> tempRecCols;
          for (const auto& recCol : recCollisions) {
            tempRecCols.push_back(std::make_pair(recCol, recCol.numContrib()));
          }
          std::sort(tempRecCols.begin(), tempRecCols.end(), compare);
          centFT0MGen = tempRecCols.at(0).first.centFT0M();
        } else {
          centFT0MGen = -999.;
        }
        registry.fill(HIST("QA/hEtaDstarGen"), mcParticle.eta());
        registry.fill(HIST("QA/hPtDstarGen"), ptGen);
        registry.fill(HIST("QA/hPtVsYDstarGen"), ptGen, yGen);
        registry.fill(HIST("Efficiency/hPtVsCentDstarGen"), ptGen, centFT0MGen);
        // only promt Dstar candidate at Generator level
        if (mcParticle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("QA/hPtPromptDstarGen"), ptGen);
          registry.fill(HIST("QA/hPtVsYPromptDstarGen"), ptGen, yGen);
        } else if (mcParticle.originMcGen() == RecoDecay::OriginType::NonPrompt) { // only non-prompt Dstar candidate at Generator level
          registry.fill(HIST("QA/hPtNonPromptDstarGen"), ptGen);
          registry.fill(HIST("QA/hPtVsYNonPromptDstarGen"), ptGen, yGen);
        }
      }
    } // MC Particle loop ends
  }

  // process data function without susing ML
  void processDataWoML(CollisionsWCent const& cols, soa::Filtered<CandDstarWSelFlag> const& selectedCands)
  {
    runTaskDstar<false, soa::Filtered<CandDstarWSelFlag>, Preslice<soa::Filtered<CandDstarWSelFlag>>>(cols, selectedCands, preslicSelectedCandDstarPerCol);
  }
  PROCESS_SWITCH(HfTaskDstarToD0Pi, processDataWoML, "Process Data without ML", true);

  // process data function with using ML, Here we store BDT score as well
  void processDataWML(CollisionsWCent const& cols, soa::Filtered<CandDstarWSelFlagWMl> const& selectedCands)
  {
    runTaskDstar<true, soa::Filtered<CandDstarWSelFlagWMl>, Preslice<soa::Filtered<CandDstarWSelFlagWMl>>>(cols, selectedCands, preslicSelectedCandDstarPerColWMl);
  }
  PROCESS_SWITCH(HfTaskDstarToD0Pi, processDataWML, "Process Data with ML", false);

  // process MC function without using ML
  void processMcWoMl(aod::McCollisions const&, CollisionsWCentMcLabel const& collisions, CandDstarWSelFlagMcRec const&,
                     CandDstarMcGen const& rowsMcPartilces,
                     aod::TracksWMc const&)
  {
    rowsSelectedCandDstarMcRec.bindExternalIndices(&collisions);
    runMcRecTaskDstar<false, Partition<CandDstarWSelFlagMcRec>>(rowsSelectedCandDstarMcRec, rowsMcPartilces);
    runMcGenTaskDstar(collisions, rowsMcPartilces);
  }
  PROCESS_SWITCH(HfTaskDstarToD0Pi, processMcWoMl, "Process MC Data without ML", false);

  // process MC function with using ML
  void processMcWML(aod::McCollisions const&, CollisionsWCentMcLabel const& collisions, CandDstarWSelFlagWMlMcRec const&,
                    CandDstarMcGen const& rowsMcPartilces,
                    aod::TracksWMc const&)
  {
    rowsSelectedCandDstarMcRecWMl.bindExternalIndices(&collisions);
    runMcRecTaskDstar<true, Partition<CandDstarWSelFlagWMlMcRec>>(rowsSelectedCandDstarMcRecWMl, rowsMcPartilces);
    runMcGenTaskDstar(collisions, rowsMcPartilces);
  }
  PROCESS_SWITCH(HfTaskDstarToD0Pi, processMcWML, "Process MC Data with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDstarToD0Pi>(cfgc)};
}
