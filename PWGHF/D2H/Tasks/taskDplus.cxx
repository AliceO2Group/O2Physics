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

/// \file taskDplus.cxx
/// \brief D± analysis task
/// \note Extended from taskD0
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, Politecnico and INFN Torino
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Luca Aglietta <luca.aglietta@cern.ch>, University and INFN Torino
/// \author Minjung Kim <minjung.kim@cern.ch>, CERN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsUpcHf.h"
#include "PWGUD/Core/UPCHelpers.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
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

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::hf_evsel;
using namespace o2::analysis::hf_upc;

/// D± analysis task
struct HfTaskDplus {
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for DPlus"}; // 7 corresponds to topo+PID cuts
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<double> yGenNBins{"yGenNBins", 100, "number of bins for y axis in sparse for gen candidates"};
  Configurable<int> centEstimator{"centEstimator", 0, "Centrality estimation (None: 0, FT0C: 2, FT0M: 3)"};
  Configurable<int> occEstimator{"occEstimator", 0, "Occupancy estimation (None: 0, ITS: 1, FT0C: 2)"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<bool> storeCentrality{"storeCentrality", false, "Flag to store centrality information"};
  Configurable<bool> storeOccupancy{"storeOccupancy", false, "Flag to store occupancy information"};
  Configurable<bool> storeIR{"storeIR", false, "Flag to store interaction rate information"};
  Configurable<bool> storePvContributors{"storePvContributors", false, "Flag to store number of PV contributors information"};
  Configurable<bool> fillMcBkgHistos{"fillMcBkgHistos", false, "Flag to fill and store histograms for MC background"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> irSource{"irSource", "ZNC hadronic", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  HfEventSelection hfEvSel;         // event selection and monitoring
  HfUpcGapThresholds upcThresholds; // UPC gap determination thresholds
  ctpRateFetcher mRateFetcher;      // interaction rate fetcher

  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandDplusDataWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandDplusMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  using CandDplusMcRecoWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec, aod::HfMlDplusToPiKPi>>;
  using CandDplusMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  using CollisionsCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs>;
  using McRecoCollisionsCent = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs>;

  Filter filterDplusFlag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi))) != static_cast<uint8_t>(0);

  Preslice<aod::HfCand3Prong> candDplusPerCollision = aod::hf_cand::collisionId;
  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<aod::McCollisionLabels> recoColPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  // data
  Partition<CandDplusData> selectedDPlusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<CandDplusDataWithMl> selectedDPlusCandidatesWithMl = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  // Matched MC
  Partition<CandDplusMcReco> recoDPlusCandidates = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<CandDplusMcRecoWithMl> recoDPlusCandidatesWithMl = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  // MC Bkg
  Partition<CandDplusMcReco> recoBkgCandidates = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<CandDplusMcRecoWithMl> recoBkgCandidatesWithMl = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {600, 1.67, 2.27}, "Cand. mass bins"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {40, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {110, 0., 110.}, "axis for centrality"};
  ConfigurableAxis thnConfigAxisOccupancy{"thnConfigAxisOccupancy", {14, 0, 14000}, "axis for occupancy"};
  ConfigurableAxis thnConfigAxisIR{"thnConfigAxisIR", {5000, 0, 500}, "Interaction rate (kHz)"};
  ConfigurableAxis thnConfigAxisPvContributors{"thnConfigAxisPvContributors", {100, 0., 100.}, "axis for PV contributors"};
  ConfigurableAxis thnConfigAxisPtBHad{"thnConfigAxisPtBHad", {25, 0., 50}, "axis for pt of B hadron decayed into D candidate"};
  ConfigurableAxis thnConfigAxisFlagBHad{"thnConfigAxisFlagBHad", {5, 0., 5}, "axis for PDG of B hadron"};
  ConfigurableAxis thnConfigAxisMlScore0{"thnConfigAxisMlScore0", {100, 0., 1.}, "axis for ML output score 0"};
  ConfigurableAxis thnConfigAxisMlScore1{"thnConfigAxisMlScore1", {100, 0., 1.}, "axis for ML output score 1"};
  ConfigurableAxis thnConfigAxisMlScore2{"thnConfigAxisMlScore2", {100, 0., 1.}, "axis for ML output score 2"};
  ConfigurableAxis thnConfigAxisGapType{"thnConfigAxisGapType", {7, -1.5, 5.5}, "axis for UPC gap type (see TrueGap enum in o2::aod::sgselector)"};
  ConfigurableAxis thnConfigAxisFT0{"thnConfigAxisFT0", {1001, -1.5, 999.5}, "axis for FT0 amplitude (a.u.)"};
  ConfigurableAxis thnConfigAxisFV0A{"thnConfigAxisFV0A", {2001, -1.5, 1999.5}, "axis for FV0-A amplitude (a.u.)"};
  ConfigurableAxis thnConfigAxisFDD{"thnConfigAxisFDD", {200, 0., 4000.}, "axis for FDD amplitude (a.u.)"};
  ConfigurableAxis thnConfigAxisZN{"thnConfigAxisZN", {510, -1.5, 49.5}, "axis for ZN energy (a.u.)"};

  HistogramRegistry registry{
    "registry",
    {{"hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hCPARecSig", "3-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecBg", "3-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hEtaRecSig", "3-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaRecBg", "3-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}}}};

  void init(InitContext&)
  {
    std::array<bool, 6> doprocess{doprocessData, doprocessDataWithMl, doprocessMc, doprocessMcWithMl, doprocessDataWithUpc, doprocessDataWithMlWithUpc};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function should be enabled! Please check your configuration!");
    }
    auto vbins = static_cast<std::vector<double>>(binsPt);
    AxisSpec const thnAxisPt = {vbins, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec const thnAxisMass = {thnConfigAxisMass, "inv. mass (K#pi#pi) (GeV/#it{c}^{2})"};
    AxisSpec const thnAxisY = {thnConfigAxisY, "y"};
    AxisSpec const thnAxisMlScore0 = {thnConfigAxisMlScore0, "Score 0"};
    AxisSpec const thnAxisMlScore1 = {thnConfigAxisMlScore1, "Score 1"};
    AxisSpec const thnAxisMlScore2 = {thnConfigAxisMlScore2, "Score 2"};
    AxisSpec const thnAxisPtBHad{thnConfigAxisPtBHad, "#it{p}_{T,B} (GeV/#it{c})"};
    AxisSpec const thnAxisFlagBHad{thnConfigAxisFlagBHad, "B Hadron flag"};
    AxisSpec const thnAxisCent{thnConfigAxisCent, "Centrality"};
    AxisSpec const thnAxisOccupancy{thnConfigAxisOccupancy, "Occupancy"};
    AxisSpec const thnAxisIR{thnConfigAxisIR, "Interaction rate (kHz)"};
    AxisSpec const thnAxisPvContributors{thnConfigAxisPvContributors, "PV contributors"};
    AxisSpec const thnAxisGapType{thnConfigAxisGapType, "Gap type"};
    AxisSpec const thnAxisFT0A{thnConfigAxisFT0, "FT0-A amplitude"};
    AxisSpec const thnAxisFT0C{thnConfigAxisFT0, "FT0-C amplitude"};
    AxisSpec const thnAxisFV0A{thnConfigAxisFV0A, "FV0-A amplitude"};
    AxisSpec const thnAxisFDDA{thnConfigAxisFDD, "FDD-A amplitude"};
    AxisSpec const thnAxisFDDC{thnConfigAxisFDD, "FDD-C amplitude"};
    AxisSpec const thnAxisZNA{thnConfigAxisZN, "ZNA energy"};
    AxisSpec const thnAxisZNC{thnConfigAxisZN, "ZNC energy"};

    registry.add("hMass", "3-prong candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{350, 1.7, 2.05}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "3-prong candidates;proper lifetime (D^{#pm}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthXY", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecayLengthXY", "3-prong candidates;norm. decay length xy;entries", {HistType::kTH2F, {{80, 0., 80.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "3-prong candidates;cos. pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxy", "3-prong candidates;cos. pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterXY", "3-prong candidates;impact parameter xy (cm);entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMaxNormalisedDeltaIP", "3-prong candidates;norm. IP;entries", {HistType::kTH2F, {{200, -20., 20.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterProngSqSum", "3-prong candidates;squared sum of prong imp. par. (cm^{2});entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthError", "3-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecayLengthXYError", "3-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpactParameterError", "3-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigPrompt", "3-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtRecBg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenSig", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtVsYRecSig_RecoPID", "3-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigPromptRecoPID", "3-prong candidates (RecoPID - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigNonPromptRecoPID", "3-prong candidates (RecoPID - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigRecoTopol", "3-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigPromptRecoTopol", "3-prong candidates (RecoTopol - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigNonPromptRecoTopol", "3-prong candidates (RecoTopol - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSig_RecoSkim", "3-prong candidates (RecoSkim - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigPrompt_RecoSkim", "3-prong candidates (RecoSkim - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYRecSigNonPrompt_RecoSkim", "3-prong candidates (RecoSkim - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYGen", "MC particles (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});
    registry.add("hPtVsYGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {100, -5., 5.}}});

    if (doprocessDataWithMl || doprocessData || doprocessDataWithMlWithUpc || doprocessDataWithUpc) {
      std::vector<AxisSpec> axes = {thnAxisMass, thnAxisPt};

      if (doprocessDataWithMl || doprocessDataWithMlWithUpc) {
        axes.push_back(thnAxisMlScore0);
        axes.push_back(thnAxisMlScore1);
        axes.push_back(thnAxisMlScore2);
      }
      if (storeCentrality) {
        axes.push_back(thnAxisCent);
      }
      if (storeOccupancy) {
        axes.push_back(thnAxisOccupancy);
      }
      if (storeIR) {
        axes.push_back(thnAxisIR);
      }
      if (doprocessDataWithMlWithUpc || doprocessDataWithUpc) {
        axes.push_back(thnAxisGapType);
        axes.push_back(thnAxisFT0A);
        axes.push_back(thnAxisFT0C);
        axes.push_back(thnAxisFV0A);
        axes.push_back(thnAxisFDDA);
        axes.push_back(thnAxisFDDC);
        axes.push_back(thnAxisZNA);
        axes.push_back(thnAxisZNC);
      }

      registry.add("hSparseMass", "THn for Dplus", HistType::kTHnSparseF, axes);
    }
    if (doprocessMcWithMl || doprocessMc) {
      std::vector<AxisSpec> axes = {thnAxisMass, thnAxisPt};
      std::vector<AxisSpec> axesFD = {thnAxisMass, thnAxisPt};
      std::vector<AxisSpec> axesGenPrompt = {thnAxisPt, thnAxisY};
      std::vector<AxisSpec> axesGenFD = {thnAxisPt, thnAxisY};

      if (doprocessMcWithMl) {
        axes.insert(axes.end(), {thnAxisMlScore0, thnAxisMlScore1, thnAxisMlScore2});
        axesFD.insert(axesFD.end(), {thnAxisMlScore0, thnAxisMlScore1, thnAxisMlScore2});
      }
      if (storeCentrality) {
        axes.push_back(thnAxisCent);
        axesFD.push_back(thnAxisCent);
        axesGenPrompt.push_back(thnAxisCent);
        axesGenFD.push_back(thnAxisCent);
      }
      if (storeOccupancy) {
        axes.push_back(thnAxisOccupancy);
        axesFD.push_back(thnAxisOccupancy);
        axesGenPrompt.push_back(thnAxisOccupancy);
        axesGenFD.push_back(thnAxisOccupancy);
      }
      if (storePvContributors) {
        axes.push_back(thnAxisPvContributors);
        axesFD.push_back(thnAxisPvContributors);
        axesGenPrompt.push_back(thnAxisPvContributors);
        axesGenFD.push_back(thnAxisPvContributors);
      }

      axesFD.push_back(thnAxisPtBHad);
      axesFD.push_back(thnAxisFlagBHad);
      axesGenFD.push_back(thnAxisPtBHad);
      axesGenFD.push_back(thnAxisFlagBHad);

      registry.add("hSparseMassPrompt", "THn for Dplus Prompt", HistType::kTHnSparseF, axes);
      registry.add("hSparseMassFD", "THn for Dplus FD", HistType::kTHnSparseF, axesFD);
      if (fillMcBkgHistos) {
        registry.add("hSparseMassBkg", "THn for Dplus Bkg", HistType::kTHnSparseF, axes);
      }
      registry.add("hSparseMassNotMatched", "THn for Dplus not matched", HistType::kTHnSparseF, axes);
      registry.add("hSparseMassGenPrompt", "THn for gen Prompt Dplus", HistType::kTHnSparseF, axesGenPrompt);
      registry.add("hSparseMassGenFD", "THn for gen FD Dplus", HistType::kTHnSparseF, axesGenFD);
    }

    registry.add("Data/fitInfo/ampFT0A_vs_ampFT0C", "FT0-A vs FT0-C amplitude;FT0-A amplitude (a.u.);FT0-C amplitude (a.u.)", {HistType::kTH2F, {{2500, 0., 250}, {2500, 0., 250}}});
    registry.add("Data/zdc/energyZNA_vs_energyZNC", "ZNA vs ZNC common energy;E_{ZNA}^{common} (a.u.);E_{ZNC}^{common} (a.u.)", {HistType::kTH2F, {{200, 0., 20}, {200, 0., 20}}});
    registry.add("Data/hUpcGapAfterSelection", "UPC gap type after selection;Gap type;Counts", {HistType::kTH1F, {{7, -1.5, 5.5}}});

    hfEvSel.addHistograms(registry); // collision monitoring

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  // Fill histograms of quantities for the reconstructed Dplus candidates
  /// \param candidate is candidate
  template <typename T1>
  void fillHisto(const T1& candidate)
  {
    float const pt = candidate.pt();
    registry.fill(HIST("hMass"), HfHelper::invMassDplusToPiKPi(candidate), pt);
    registry.fill(HIST("hPt"), pt);
    registry.fill(HIST("hEta"), candidate.eta(), pt);
    registry.fill(HIST("hCt"), HfHelper::ctDplus(candidate), pt);
    registry.fill(HIST("hDecayLength"), candidate.decayLength(), pt);
    registry.fill(HIST("hDecayLengthXY"), candidate.decayLengthXY(), pt);
    registry.fill(HIST("hNormalisedDecayLengthXY"), candidate.decayLengthXYNormalised(), pt);
    registry.fill(HIST("hCPA"), candidate.cpa(), pt);
    registry.fill(HIST("hCPAxy"), candidate.cpaXY(), pt);
    registry.fill(HIST("hImpactParameterXY"), candidate.impactParameterXY(), pt);
    registry.fill(HIST("hMaxNormalisedDeltaIP"), candidate.maxNormalisedDeltaIP(), pt);
    registry.fill(HIST("hImpactParameterProngSqSum"), candidate.impactParameterProngSqSum(), pt);
    registry.fill(HIST("hDecayLengthError"), candidate.errorDecayLength(), pt);
    registry.fill(HIST("hDecayLengthXYError"), candidate.errorDecayLengthXY(), pt);
    registry.fill(HIST("hImpactParameterError"), candidate.errorImpactParameter0(), pt);
    registry.fill(HIST("hImpactParameterError"), candidate.errorImpactParameter1(), pt);
    registry.fill(HIST("hImpactParameterError"), candidate.errorImpactParameter2(), pt);
    registry.fill(HIST("hPtProng0"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1"), candidate.ptProng1());
    registry.fill(HIST("hPtProng2"), candidate.ptProng2());
    registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), pt);
    registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), pt);
    registry.fill(HIST("hd0Prong2"), candidate.impactParameter2(), pt);
  }

  // Fill THnSparses for the ML analysis
  /// \param candidate is a particle candidate
  /// \param ptbhad transverse momentum of beauty mother for nonprompt candidates
  /// \param flagBHad transverse momentum of beauty mother for nonprompt candidates
  /// \param centrality collision centrality
  /// \param occupancy collision occupancy
  /// \param numPvContributors contributors to the PV
  template <bool IsMc, bool IsMatched, typename T1>
  void fillSparseML(const T1& candidate,
                    float ptbhad,
                    int flagBHad,
                    float centrality,
                    float occupancy,
                    float numPvContributors)
  {
    std::vector<float> outputMl = {-999., -999., -999.};
    for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
      outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
    }
    if constexpr (IsMc) {                                               // MC
      if constexpr (IsMatched) {                                        // Matched
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) { // Prompt

          if (storeCentrality && storeOccupancy) {
            registry.fill(HIST("hSparseMassPrompt"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], centrality, occupancy);
          } else if (storeCentrality && !storeOccupancy) {
            registry.fill(HIST("hSparseMassPrompt"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], centrality);
          } else if (!storeCentrality && storeOccupancy) {
            registry.fill(HIST("hSparseMassPrompt"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], occupancy);
          } else if (!storeCentrality && !storeOccupancy && storePvContributors) {
            registry.fill(HIST("hSparseMassPrompt"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], numPvContributors);
          } else {
            registry.fill(HIST("hSparseMassPrompt"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2]);
          }

        } else if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) { // FD

          if (storeCentrality && storeOccupancy) {
            registry.fill(HIST("hSparseMassFD"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], centrality, occupancy, ptbhad, flagBHad);
          } else if (storeCentrality && !storeOccupancy) {
            registry.fill(HIST("hSparseMassFD"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], centrality, ptbhad, flagBHad);
          } else if (!storeCentrality && storeOccupancy) {
            registry.fill(HIST("hSparseMassFD"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], occupancy, ptbhad, flagBHad);
          } else if (!storeCentrality && !storeOccupancy && storePvContributors) {
            registry.fill(HIST("hSparseMassFD"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], numPvContributors, ptbhad, flagBHad);
          } else {
            registry.fill(HIST("hSparseMassFD"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], ptbhad, flagBHad);
          }

        } else { // Bkg
          if (fillMcBkgHistos) {
            if (storeCentrality && storeOccupancy) {
              registry.fill(HIST("hSparseMassBkg"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], centrality, occupancy);
            } else if (storeCentrality && !storeOccupancy) {
              registry.fill(HIST("hSparseMassBkg"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], centrality);
            } else if (!storeCentrality && storeOccupancy) {
              registry.fill(HIST("hSparseMassBkg"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], occupancy);
            } else if (!storeCentrality && !storeOccupancy && storePvContributors) {
              registry.fill(HIST("hSparseMassBkg"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], numPvContributors);
            } else {
              registry.fill(HIST("hSparseMassBkg"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2]);
            }
          }
        }
      } else {
        if (storeCentrality && storeOccupancy) {
          registry.fill(HIST("hSparseMassNotMatched"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], centrality, occupancy);
        } else if (storeCentrality && !storeOccupancy) {
          registry.fill(HIST("hSparseMassNotMatched"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], centrality);
        } else if (!storeCentrality && storeOccupancy) {
          registry.fill(HIST("hSparseMassNotMatched"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], occupancy);
        } else if (!storeCentrality && !storeOccupancy && storePvContributors) {
          registry.fill(HIST("hSparseMassNotMatched"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], numPvContributors);
        } else {
          registry.fill(HIST("hSparseMassNotMatched"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2]);
        }
      }
    } else { // Data
      if (storeCentrality && storeOccupancy) {
        registry.fill(HIST("hSparseMass"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], centrality, occupancy);
      } else if (storeCentrality && !storeOccupancy) {
        registry.fill(HIST("hSparseMass"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], centrality);
      } else if (!storeCentrality && storeOccupancy) {
        registry.fill(HIST("hSparseMass"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], occupancy);
      } else if (!storeCentrality && !storeOccupancy && storePvContributors) {
        registry.fill(HIST("hSparseMass"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2], numPvContributors);
      } else {
        registry.fill(HIST("hSparseMass"), HfHelper::invMassDplusToPiKPi(candidate), candidate.pt(), outputMl[0], outputMl[1], outputMl[2]);
      }
    }
  }

  // Fill histograms of quantities for the reconstructed Dplus candidates with MC matching
  /// \param candidate is candidate
  template <bool IsMatched, typename T1>
  void fillHistoMCRec(const T1& candidate)
  {
    if constexpr (IsMatched) {
      auto ptRec = candidate.pt();
      auto yRec = HfHelper::yDplus(candidate);
      registry.fill(HIST("hPtVsYRecSig_RecoSkim"), ptRec, yRec);
      if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoTopol)) {
        registry.fill(HIST("hPtVsYRecSigRecoTopol"), ptRec, yRec);
      }
      if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoPID)) {
        registry.fill(HIST("hPtVsYRecSig_RecoPID"), ptRec, yRec);
      }
      if (candidate.isSelDplusToPiKPi() >= selectionFlagDplus) {
        registry.fill(HIST("hPtRecSig"), ptRec); // rec. level pT
      }
      if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hPtVsYRecSigPrompt_RecoSkim"), ptRec, yRec);
        if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoTopol)) {
          registry.fill(HIST("hPtVsYRecSigPromptRecoTopol"), ptRec, yRec);
        }
        if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoPID)) {
          registry.fill(HIST("hPtVsYRecSigPromptRecoPID"), ptRec, yRec);
        }
        if (candidate.isSelDplusToPiKPi() >= selectionFlagDplus) {
          registry.fill(HIST("hPtRecSigPrompt"), ptRec); // rec. level pT, prompt
        }
      } else {
        registry.fill(HIST("hPtVsYRecSigNonPrompt_RecoSkim"), ptRec, yRec);
        if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoTopol)) {
          registry.fill(HIST("hPtVsYRecSigNonPromptRecoTopol"), ptRec, yRec);
        }
        if (TESTBIT(candidate.isSelDplusToPiKPi(), aod::SelectionStep::RecoPID)) {
          registry.fill(HIST("hPtVsYRecSigNonPromptRecoPID"), ptRec, yRec);
        }
        if (candidate.isSelDplusToPiKPi() >= selectionFlagDplus) {
          registry.fill(HIST("hPtRecSigNonPrompt"), ptRec); // rec. level pT, non-prompt
        }
      }
      registry.fill(HIST("hCPARecSig"), candidate.cpa());
      registry.fill(HIST("hEtaRecSig"), candidate.eta());
    } else {
      registry.fill(HIST("hPtRecBg"), candidate.pt());
      registry.fill(HIST("hCPARecBg"), candidate.cpa());
      registry.fill(HIST("hEtaRecBg"), candidate.eta());
    }
  }

  // Fill histograms of quantities for generated Dplus particles
  /// \param particle is a particle with MC information
  template <typename T2>
  void fillHistoMCGen(const T2& particle)
  {
    auto ptGen = particle.pt();
    auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDPlus);
    registry.fill(HIST("hPtGen"), ptGen);
    registry.fill(HIST("hPtVsYGen"), ptGen, yGen);
    if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtGenPrompt"), ptGen);
      registry.fill(HIST("hPtVsYGenPrompt"), ptGen, yGen);
    } else {
      registry.fill(HIST("hPtGenNonPrompt"), ptGen);
      registry.fill(HIST("hPtVsYGenNonPrompt"), ptGen, yGen);
    }
    registry.fill(HIST("hEtaGen"), particle.eta());
  }

  // Fill THnSparse of quantities for generated Dplus particles
  /// \param particle is a particle with MC information
  /// \param ptGenB transverse momentum of beauty mother for nonprompt candidates
  /// \param flagGenB transverse momentum of beauty mother for nonprompt candidates
  /// \param centrality collision centrality
  /// \param occupancy collision occupancy
  /// \param numPvContributors contributors to the PV
  template <typename T1>
  void fillSparseMcGen(const T1& particle,
                       float ptGenB,
                       int flagGenB,
                       float centrality,
                       float occupancy,
                       float numPvContributors)
  {
    auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDPlus);
    if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
      if (storeCentrality && storeOccupancy) {
        registry.fill(HIST("hSparseMassGenPrompt"), particle.pt(), yGen, centrality, occupancy);
      } else if (storeCentrality && !storeOccupancy) {
        registry.fill(HIST("hSparseMassGenPrompt"), particle.pt(), yGen, centrality);
      } else if (!storeCentrality && storeOccupancy) {
        registry.fill(HIST("hSparseMassGenPrompt"), particle.pt(), yGen, occupancy);
      } else if (!storeCentrality && !storeOccupancy && storePvContributors) {
        registry.fill(HIST("hSparseMassGenPrompt"), particle.pt(), yGen, numPvContributors);
      } else {
        registry.fill(HIST("hSparseMassGenPrompt"), particle.pt(), yGen);
      }
    } else {
      if (storeCentrality && storeOccupancy) {
        registry.fill(HIST("hSparseMassGenFD"), particle.pt(), yGen, centrality, occupancy, ptGenB, flagGenB);
      } else if (storeCentrality && !storeOccupancy) {
        registry.fill(HIST("hSparseMassGenFD"), particle.pt(), yGen, centrality, ptGenB, flagGenB);
      } else if (!storeCentrality && storeOccupancy) {
        registry.fill(HIST("hSparseMassGenFD"), particle.pt(), yGen, occupancy, ptGenB, flagGenB);
      } else if (!storeCentrality && !storeOccupancy && storePvContributors) {
        registry.fill(HIST("hSparseMassGenFD"), particle.pt(), yGen, numPvContributors, ptGenB, flagGenB);
      } else {
        registry.fill(HIST("hSparseMassGenFD"), particle.pt(), yGen, ptGenB, flagGenB);
      }
    }
  }

  // Run analysis for the reconstructed Dplus candidates from data
  /// \param candidates are reconstructed candidates
  template <bool FillMl, typename T1>
  void runDataAnalysis(const T1& /*candidates*/, CollisionsCent const& /*colls*/)
  {
    float cent{-1.f};
    float occ{-1.f};
    float numPvContr{-1.f};
    float ptBhad{-1.f};
    int const flagBHad{-1};
    if constexpr (!FillMl) {
      for (const auto& candidate : selectedDPlusCandidates) {
        if ((yCandRecoMax >= 0. && std::abs(HfHelper::yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }
        fillHisto(candidate);
      }
    } else {
      for (const auto& candidate : selectedDPlusCandidatesWithMl) {
        if ((yCandRecoMax >= 0. && std::abs(HfHelper::yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }

        if (storeCentrality || storeOccupancy) {
          auto collision = candidate.template collision_as<CollisionsCent>();
          if (storeCentrality && centEstimator != CentralityEstimator::None) {
            cent = getCentralityColl(collision, centEstimator);
          }
          if (storeOccupancy && occEstimator != OccupancyEstimator::None) {
            occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
          }
          if (storePvContributors) {
            numPvContr = collision.numContrib();
          }
        }

        fillHisto(candidate);
        fillSparseML<false, false>(candidate, ptBhad, flagBHad, cent, occ, numPvContr);
      }
    }
  }

  // Run analysis for the reconstructed Dplus candidates with MC matching
  /// \param recoCandidates are reconstructed candidates
  /// \param recoColls are reconstructed collisions
  template <bool FillMl>
  void runAnalysisMcRec(McRecoCollisionsCent const& /*recoColls*/)
  {
    float cent{-1};
    float occ{-1};
    float numPvContr{-1};
    float ptBhad{-1};
    int flagBHad{-1};

    // MC rec. w/o Ml
    if constexpr (!FillMl) {
      for (const auto& candidate : recoDPlusCandidates) {
        if ((yCandRecoMax >= 0. && std::abs(HfHelper::yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }
        fillHisto(candidate);
        fillHistoMCRec<true>(candidate);
      }
      // Bkg
      if (fillMcBkgHistos) {
        for (const auto& candidate : recoBkgCandidates) {
          if ((yCandRecoMax >= 0. && std::abs(HfHelper::yDplus(candidate)) > yCandRecoMax)) {
            continue;
          }
          fillHistoMCRec<false>(candidate);
        }
      }
    } else {
      for (const auto& candidate : recoDPlusCandidatesWithMl) {
        if ((yCandRecoMax >= 0. && std::abs(HfHelper::yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }
        ptBhad = candidate.ptBhadMotherPart();
        flagBHad = getBHadMotherFlag(candidate.pdgBhadMotherPart());
        auto collision = candidate.template collision_as<McRecoCollisionsCent>();
        if (storeCentrality || storeOccupancy) {
          if (storeCentrality && centEstimator != CentralityEstimator::None) {
            cent = getCentralityColl(collision, centEstimator);
          }
          if (storeOccupancy && occEstimator != OccupancyEstimator::None) {
            occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
          }
        }
        if (storePvContributors) {
          numPvContr = collision.numContrib();
        }
        fillHisto(candidate);
        fillHistoMCRec<true>(candidate);
        fillSparseML<true, true>(candidate, ptBhad, flagBHad, cent, occ, numPvContr);
      }
      // Bkg
      ptBhad = -1;
      flagBHad = -1;
      if (fillMcBkgHistos) {
        for (const auto& candidate : recoBkgCandidatesWithMl) {
          if ((yCandRecoMax >= 0. && std::abs(HfHelper::yDplus(candidate)) > yCandRecoMax)) {
            continue;
          }
          auto collision = candidate.template collision_as<McRecoCollisionsCent>();
          if (storeCentrality && centEstimator != CentralityEstimator::None) {
            cent = getCentralityColl(collision, centEstimator);
          }
          if (storeOccupancy && occEstimator != OccupancyEstimator::None) {
            occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
          }
          if (storePvContributors) {
            numPvContr = collision.numContrib();
          }
          fillHistoMCRec<false>(candidate);
          fillSparseML<true, false>(candidate, ptBhad, flagBHad, cent, occ, numPvContr);
        }
      }
    }
  }

  // Run analysis for the generated Dplus candidates
  /// \param mcGenCollisions are the generated MC collisions
  /// \param mcRecoCollisions are the reconstructed MC collisions
  /// \param mcGenParticles are the generated MC particle candidates
  template <bool FillMl, typename Cand>
  void runAnalysisMcGen(aod::McCollisions const& mcGenCollisions,
                        McRecoCollisionsCent const& mcRecoCollisions,
                        Cand const& mcGenParticles)
  {
    // MC gen.
    float cent{-1.};
    float occ{-1.};
    float numPvContr{-1.};
    float ptGenB{-1.};
    int flagGenB{-1};

    for (const auto& mcGenCollision : mcGenCollisions) {
      const auto recoCollsPerGenMcColl = mcRecoCollisions.sliceBy(recoColPerMcCollision, mcGenCollision.globalIndex());
      const auto mcParticlesPerGenMcColl = mcGenParticles.sliceBy(mcParticlesPerMcCollision, mcGenCollision.globalIndex());
      if (storeCentrality && centEstimator != CentralityEstimator::None) {
        cent = getCentralityGenColl(recoCollsPerGenMcColl, centEstimator);
      }
      if (storeOccupancy && occEstimator != OccupancyEstimator::None) {
        occ = o2::hf_occupancy::getOccupancyGenColl(recoCollsPerGenMcColl, occEstimator);
      }

      for (const auto& particle : mcParticlesPerGenMcColl) {
        ptGenB = -1;
        flagGenB = -1;
        numPvContr = -1;
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassDPlus);
        if ((yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) || (std::abs(particle.flagMcMatchGen()) != hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)) {
          continue;
        }
        if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          auto bHadMother = mcGenParticles.rawIteratorAt(particle.idxBhadMotherPart() - mcGenParticles.offset());
          flagGenB = getBHadMotherFlag(bHadMother.pdgCode());
          ptGenB = bHadMother.pt();
        }
        for (const auto& recCol : recoCollsPerGenMcColl) {
          numPvContr = std::max<float>(numPvContr, recCol.numContrib());
        }
        fillHistoMCGen(particle);
        if constexpr (FillMl) {
          fillSparseMcGen(particle, ptGenB, flagGenB, cent, occ, numPvContr);
        }
      }
    }
  }

  template <bool FillMl, typename CollType, typename CandType, typename BCsType>
  void runAnalysisPerCollisionDataWithUpc(CollType const& collisions,
                                          CandType const& candidates,
                                          BCsType const& bcs,
                                          aod::FT0s const& ft0s,
                                          aod::FV0As const& fv0as,
                                          aod::FDDs const& fdds)
  {
    for (const auto& collision : collisions) {
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMaskWithUpc<true, CentralityEstimator::None, BCsType>(collision, centrality, ccdb, registry, bcs);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }
      const auto& bc = collision.template bc_as<BCsType>();

      // Determine gap type using SGSelector with BC range checking
      const auto gapResult = hf_upc::determineGapType(collision, bcs, upcThresholds);
      const int gap = gapResult.value;

      // Use the BC with FIT activity if available from SGSelector
      auto bcForUPC = bc;
      if (gapResult.bc) {
        bcForUPC = *(gapResult.bc);
      }

      // Get FIT information from the UPC BC
      upchelpers::FITInfo fitInfo{};
      udhelpers::getFITinfo(fitInfo, bcForUPC, bcs, ft0s, fv0as, fdds);

      // Get ZDC energies if available (extract once and reuse)
      const bool hasZdc = bcForUPC.has_zdc();
      float zdcEnergyZNA = -1.f;
      float zdcEnergyZNC = -1.f;
      if (hasZdc) {
        const auto& zdc = bcForUPC.zdc();
        zdcEnergyZNA = zdc.energyCommonZNA();
        zdcEnergyZNC = zdc.energyCommonZNC();
        registry.fill(HIST("Data/zdc/energyZNA_vs_energyZNC"), zdcEnergyZNA, zdcEnergyZNC);
      }
      registry.fill(HIST("Data/fitInfo/ampFT0A_vs_ampFT0C"), fitInfo.ampFT0A, fitInfo.ampFT0C);
      registry.fill(HIST("Data/hUpcGapAfterSelection"), gap);

      const auto thisCollId = collision.globalIndex();
      const auto& groupedDplusCandidates = candidates.sliceBy(candDplusPerCollision, thisCollId);

      float cent{-1.f};
      float occ{-1.f};
      float ir{-1.f};
      if (storeOccupancy && occEstimator != OccupancyEstimator::None) {
        occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
      }
      if (storeIR) {
        ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource, true) * 1.e-3; // kHz
      }
      static constexpr auto HSparseMass = HIST("hSparseMass");
      // Lambda function to fill THn - handles both ML and non-ML cases
      auto fillTHnData = [&](const auto& candidate) {
        // Pre-calculate vector size to avoid reallocations
        constexpr int NAxesBase = 10;                  // mass, pt, gapType, FT0A, FT0C, FV0A, FDDA, FDDC, ZNA, ZNC
        constexpr int NAxesMl = FillMl ? 3 : 0;        // 3 ML scores if FillMl
        int const nAxesCent = storeCentrality ? 1 : 0; // centrality if storeCentrality
        int const nAxesOcc = storeOccupancy ? 1 : 0;   // occupancy if storeOccupancy
        int const nAxesIR = storeIR ? 1 : 0;           // IR if storeIR
        int const nAxesTotal = NAxesBase + NAxesMl + nAxesCent + nAxesOcc + nAxesIR;

        std::vector<double> valuesToFill;
        valuesToFill.reserve(nAxesTotal);

        // Fill values in order matching histogram axes
        valuesToFill.push_back(HfHelper::invMassDplusToPiKPi(candidate));
        valuesToFill.push_back(static_cast<double>(candidate.pt()));
        if constexpr (FillMl) {
          std::vector<float> outputMl = {-999., -999., -999.};
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
          }
          valuesToFill.push_back(outputMl[0]);
          valuesToFill.push_back(outputMl[1]);
          valuesToFill.push_back(outputMl[2]);
        }
        if (storeCentrality) {
          valuesToFill.push_back(cent);
        }
        if (storeOccupancy) {
          valuesToFill.push_back(occ);
        }
        if (storeIR) {
          valuesToFill.push_back(ir);
        }
        valuesToFill.push_back(static_cast<double>(gap));
        valuesToFill.push_back(static_cast<double>(fitInfo.ampFT0A));
        valuesToFill.push_back(static_cast<double>(fitInfo.ampFT0C));
        valuesToFill.push_back(static_cast<double>(fitInfo.ampFV0A));
        valuesToFill.push_back(static_cast<double>(fitInfo.ampFDDA));
        valuesToFill.push_back(static_cast<double>(fitInfo.ampFDDC));
        valuesToFill.push_back(static_cast<double>(zdcEnergyZNA));
        valuesToFill.push_back(static_cast<double>(zdcEnergyZNC));
        registry.get<THnSparse>(HSparseMass)->Fill(valuesToFill.data());
      };

      for (const auto& candidate : groupedDplusCandidates) {
        if ((yCandRecoMax >= 0. && std::abs(HfHelper::yDplus(candidate)) > yCandRecoMax)) {
          continue;
        }
        fillHisto(candidate);
        fillTHnData(candidate);
      }
    }
  }

  // process functions
  void processData(CandDplusData const& candidates, CollisionsCent const& collisions)
  {
    runDataAnalysis<false>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskDplus, processData, "Process data w/o ML", true);

  void processDataWithMl(CandDplusDataWithMl const& candidates, CollisionsCent const& collisions)
  {
    runDataAnalysis<true>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskDplus, processDataWithMl, "Process data with ML", false);

  void processMc(CandDplusMcReco const&,
                 CandDplusMcGen const& mcGenParticles,
                 McRecoCollisionsCent const& mcRecoCollisions,
                 aod::McCollisions const& mcGenCollisions)
  {
    runAnalysisMcRec<false>(mcRecoCollisions);
    runAnalysisMcGen<false>(mcGenCollisions, mcRecoCollisions, mcGenParticles);
  }
  PROCESS_SWITCH(HfTaskDplus, processMc, "Process MC w/o ML", false);

  void processMcWithMl(CandDplusMcRecoWithMl const&,
                       CandDplusMcGen const& mcGenParticles,
                       McRecoCollisionsCent const& mcRecoCollisions,
                       aod::McCollisions const& mcGenCollisions)
  {
    runAnalysisMcRec<true>(mcRecoCollisions);
    runAnalysisMcGen<true>(mcGenCollisions, mcRecoCollisions, mcGenParticles);
  }
  PROCESS_SWITCH(HfTaskDplus, processMcWithMl, "Process MC with ML", false);

  void processDataWithUpc(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                          aod::BcFullInfos const& bcs,
                          CandDplusData const& selectedDplusCandidates,
                          aod::Tracks const&,
                          aod::FT0s const& ft0s,
                          aod::FV0As const& fv0as,
                          aod::FDDs const& fdds,
                          aod::Zdcs const& /*zdcs*/)
  {
    runAnalysisPerCollisionDataWithUpc<false>(collisions, selectedDplusCandidates, bcs, ft0s, fv0as, fdds);
  }
  PROCESS_SWITCH(HfTaskDplus, processDataWithUpc, "Process real data w/o ML with UPC", false);

  void processDataWithMlWithUpc(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                aod::BcFullInfos const& bcs,
                                CandDplusDataWithMl const& selectedDplusCandidatesMl,
                                aod::Tracks const&,
                                aod::FT0s const& ft0s,
                                aod::FV0As const& fv0as,
                                aod::FDDs const& fdds,
                                aod::Zdcs const& /*zdcs*/)
  {
    runAnalysisPerCollisionDataWithUpc<true>(collisions, selectedDplusCandidatesMl, bcs, ft0s, fv0as, fdds);
  }
  PROCESS_SWITCH(HfTaskDplus, processDataWithMlWithUpc, "Process real data with the ML method with UPC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDplus>(cfgc)};
}
