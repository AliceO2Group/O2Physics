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
///
/// \brief Analysis of D0/Lambda_c yield as a function of flattenicity
/// \author Laszlo Gyulai, laszlo.gyulai@cern.ch

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;

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
}

static const int nCellsFV0 = 48;
static const int CinnerFV0 = 32;
std::array<float, nCellsFV0> rhoLatticeFV0{0};
std::array<float, nCellsFV0> fv0AmplitudeWoCalib{0};
float calib[48] = {1.01697, 1.122, 1.03854, 1.108, 1.11634, 1.14971, 1.19321, 1.06866, 0.954675, 0.952695, 0.969853, 0.957557, 0.989784, 1.01549, 1.02182, 0.976005, 1.01865, 1.06871, 1.06264, 1.02969, 1.07378, 1.06622, 1.15057, 1.0433, 0.83654, 0.847178, 0.890027, 0.920814, 0.888271, 1.04662, 0.8869, 0.856348, 0.863181, 0.906312, 0.902166, 1.00122, 1.03303, 0.887866, 0.892437, 0.906278, 0.884976, 0.864251, 0.917221, 1.10618, 1.04028, 0.893184, 0.915734, 0.892676};

struct FlattenicityDLc {
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<int> selectionTopol{"selectionTopol", 1, "Selection Flag for topologically selected candidates"};
  Configurable<int> selectionCand{"selectionCand", 1, "Selection Flag for conj. topol. selected candidates"};
  Configurable<int> selectionPid{"selectionPid", 1, "Selection Flag for reco PID candidates"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int> selectionFlagHf{"selectionFlagHf", 1, "Selection Flag for HF flagged candidates"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};

  HfEventSelection hfEvSel;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  using TracksWPid = soa::Join<o2::aod::FullTracks, aod::TracksDCA, o2::aod::TrackSelection, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using TracksSelQuality = soa::Join<aod::TracksExtra, aod::TracksWMc>;

  using D0Candidates = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  using D0CandidatesMc = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>;

  using LcCandidates = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using LcCandidatesMc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>;

  Filter filterD0 = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK))) != static_cast<uint8_t>(0);
  Filter filterLc = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc;

  Preslice<aod::HfCand2Prong> candD0PerCollision = aod::hf_cand::collisionId;
  Preslice<aod::HfCand3Prong> candLcPerCollision = aod::hf_cand::collisionId;
  PresliceUnsorted<CollisionsWithMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollisionLc = aod::mcparticle::mcCollisionId;

  Partition<D0Candidates> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesMc> selectedD0CandidatesMc = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;

  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {120, 1.5848, 2.1848}, "Cand. inv-mass bins"};
  ConfigurableAxis thnConfigAxisFlat{"thnConfigAxisFlat", {100, 0, 1}, "Flattenicity in event associated with candidate"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {500, 0, 50}, "Cand. pT bins"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisCandType{"thnConfigAxisCandType", {6, -0.5, 5.5}, "D0 type"};
  ConfigurableAxis thnConfigAxisPtB{"thnConfigAxisPtB", {1000, 0, 100}, "Cand. beauty mother pTB bins"};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, "Cand. origin type"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  ConfigurableAxis thnConfigAxisGenPtD{"thnConfigAxisGenPtD", {500, 0, 50}, "Gen Pt D"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {1000, 0, 100}, "Gen Pt B"};
  ConfigurableAxis thnConfigAxisMassLc{"thnConfigAxisMassLc", {300, 1.98, 2.58}, ""};
  ConfigurableAxis thnConfigAxisCanType{"thnConfigAxisCanType", {5, 0., 5.}, ""};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(InitContext const&)
  {
    std::array<bool, 3> doprocess{doprocessData, doprocessMCD0, doprocessMCLc};

    hfEvSel.addHistograms(registry);
    registry.add("Flattenicity", "Number of events; 1-#rho; Counter", {kTH1F, {{100, 0, 1}}});
    registry.add("Flattenicity_calibrated", "Number of events; 1-#rho; Counter", {kTH1F, {{100, 0, 1}}});

    auto vbins = (std::vector<double>)binsPt;

    const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (#pi K) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisFlat{thnConfigAxisFlat, "#it{#rho}"};
    const AxisSpec thnAxisY{thnConfigAxisY, "y"};
    const AxisSpec thnAxisCandType{thnConfigAxisCandType, "D0 type"};
    const AxisSpec thnAxisPtB{thnConfigAxisPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisOrigin{thnConfigAxisOrigin, "Origin"};
    const AxisSpec thnAxisGenPtD{thnConfigAxisGenPtD, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisGenPtB{thnConfigAxisGenPtB, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisMassLc{thnConfigAxisMassLc, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisCanType{thnConfigAxisCanType, "candidates type"};

    std::vector<AxisSpec> axesD0 = {thnAxisMass, thnAxisPt, thnAxisFlat, thnAxisY, thnAxisCandType};
    std::vector<AxisSpec> axesD0Gen = {thnAxisGenPtD, thnAxisFlat, thnAxisGenPtB, thnAxisY, thnAxisOrigin};
    std::vector<AxisSpec> axesLcData = {thnAxisMassLc, thnAxisPt, thnAxisFlat};
    std::vector<AxisSpec> axesLcMc = {thnAxisMassLc, thnAxisPt, thnAxisFlat, thnAxisPtB, thnAxisCanType};
    std::vector<AxisSpec> axesGenLc = {thnAxisPt, thnAxisFlat, thnAxisPtB, thnAxisCanType};

    if (doprocessData) {
      registry.add("Data/D0/hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hMassFinerBinning", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hMassVsPhi", "2-prong candidates vs phi;inv. mass (#pi K) (GeV/#it{c}^{2});phi (rad);entries", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {32, 0, o2::constants::math::TwoPI}}});
      registry.add("Data/D0/hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("Data/D0/hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("Data/D0/hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("Data/D0/hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hNormalisedDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hNormalisedDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hd0ErrProng0", "2-prong candidates;prong 0 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hd0ErrProng1", "2-prong candidates;prong 1 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hDecLengthFinerBinning", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{400, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hDecLengthxyFinerBinning", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{400, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hd0Prong0FinerBinning", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hd0Prong1FinerBinning", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hd0d0FinerBinning", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -0.1, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hCTSFinerBinning", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hCtFinerBinning", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{500, -0., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hCPAFinerBinning", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/D0/hCPAXYFinerBinning", "2-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

      registry.add("Data/D0/hMassVsPtVsFlatVsYVsD0Type", "Thn for D0 candidates", HistType::kTHnSparseD, axesD0);

      registry.add("Data/Lc/hMass", ";inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{600, 1.98, 2.58}}});
      registry.add("Data/Lc/hMassVsPt", "Inv.mass vs pt;inv. mass (p K #pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins}}});
      registry.add("Data/Lc/hPt", "#pT of candidates;it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("Data/Lc/hPtProng0", "pT of prong 0;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("Data/Lc/hPtProng1", "pT of prong 1;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("Data/Lc/hPtProng2", "pT of prong 2;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("Data/Lc/hd0Prong0", ";prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("Data/Lc/hd0Prong1", ";prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("Data/Lc/hd0Prong2", ";prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("Data/Lc/hd0VsPtProng0", ";prong 0 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("Data/Lc/hd0VsPtProng1", ";prong 1 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("Data/Lc/hd0VsPtProng2", ";prong 2 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("Data/Lc/hDecLength", ";decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}});
      registry.add("Data/Lc/hDecLengthVsPt", ";decay length (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});
      registry.add("Data/Lc/hDecLengthxy", ";decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}});
      registry.add("Data/Lc/hDecLengthxyVsPt", ";decay length xy (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});
      registry.add("Data/Lc/hCt", ";proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}});
      registry.add("Data/Lc/hCtVsPt", ";proper lifetime (#Lambda_{c}) * #it{c} (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 0.2}, {vbins}}});
      registry.add("Data/Lc/hCPA", ";cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("Data/Lc/hCPAVsPt", ";cosine of pointing angle;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});
      registry.add("Data/Lc/hCPAxy", ";cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("Data/Lc/hCPAxyVsPt", ";cosine of pointing angle xy;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});
      registry.add("Data/Lc/hEta", ";#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("Data/Lc/hEtaVsPt", ";candidate #it{#eta};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
      registry.add("Data/Lc/hPhi", "#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}});
      registry.add("Data/Lc/hPhiVsPt", ";candidate #it{#Phi};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});
      registry.add("Data/Lc/hSelectionStatus", "3-prong candidates;selection status;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("Data/Lc/hImpParErrProng0VsPt", ";prong 0 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("Data/Lc/hImpParErrProng1VsPt", ";prong 1 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("Data/Lc/hImpParErrProng2VsPt", ";prong 2 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("Data/Lc/hDecLenErrVsPt", ";decay length error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 1.}, {vbins}}});

      registry.add("Data/Lc/hLcVsPtVsFlat", "THn for Reconstructed Lambdac candidates for data", HistType::kTHnSparseF, axesLcData);
    }

    if (doprocessMCD0) {
      registry.add("MC/D0/hPtGenSig", "2-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/D0/hPtVsYRecSigRecoHFFlag", "2-prong candidates (RecoHFFlag - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtGenVsPtRecSig", "2-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH2F, {{360, 0., 36.}, {360, 0., 36.}}});
      registry.add("MC/D0/hYGenVsYRecSig", "2-prong candidates (matched);#it{y}^{gen.} ;#it{y}^{rec.} ;entries", {HistType::kTH2F, {{300, -1.5, 1.5}, {300, -1.5, 1.5}}});
      registry.add("MC/D0/hMassVsPtGenVsPtRecSig", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{gen.} (GeV/#it{c});#it{p}_{T}^{rec.} (GeV/#it{c})", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {150, 0., 30.}}});
      registry.add("MC/D0/hCPARecSig", "2-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("MC/D0/hEtaRecSig", "2-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}});
      registry.add("MC/D0/hPtRecBg", "2-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/D0/hCPARecBg", "2-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("MC/D0/hEtaRecBg", "2-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}});
      registry.add("MC/D0/hMassSigBkgD0", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}});
      registry.add("MC/D0/hPtProng0Sig", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}});
      registry.add("MC/D0/hPtProng1Sig", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}});
      registry.add("MC/D0/hDecLengthSig", "2-prong candidates (matched);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}});
      registry.add("MC/D0/hDecLengthXYSig", "2-prong candidates (matched);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}});
      registry.add("MC/D0/hNormalisedDecLengthSig", "2-prong candidates (matched);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}});
      registry.add("MC/D0/hNormalisedDecLengthXYSig", "2-prong candidates (matched);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}});
      registry.add("MC/D0/hd0Prong0Sig", "2-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}});
      registry.add("MC/D0/hd0Prong1Sig", "2-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}});
      registry.add("MC/D0/hd0d0Sig", "2-prong candidates (matched);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}});
      registry.add("MC/D0/hCTSSig", "2-prong candidates (matched);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}});
      registry.add("MC/D0/hCtSig", "2-prong candidates (matched);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}});
      registry.add("MC/D0/hCPASig", "2-prong candidates (matched);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}});
      registry.add("MC/D0/hCPAxySig", "2-prong candidates (matched);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}});
      registry.add("MC/D0/hd0Prong0VsPtSig", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/D0/hd0Prong1VsPtSig", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/D0/hd0d0VsPtSig", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2}) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -0.1, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/D0/hCTSVsPtSig", "2-prong candidates;cos #it{#theta}* (D^{0}) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/D0/hCPAVsPtSig", "2-prong candidates;cosine of pointing angle vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/D0/hCPAXYVsPtSig", "2-prong candidates;cosine of pointing angle xy vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/D0/hNormalisedDecLengthxyVsPtSig", "2-prong candidates;decay length xy (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/D0/hDecLengthVsPtSig", "2-prong candidates;decay length (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/D0/hDecLengthxyVsPtSig", "2-prong candidates;decay length xy (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("MC/D0/hMassSigD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}});
      registry.add("MC/D0/hPtProng0Bkg", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}});
      registry.add("MC/D0/hPtProng1Bkg", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}});
      registry.add("MC/D0/hDecLengthBkg", "2-prong candidates (checked);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}});
      registry.add("MC/D0/hDecLengthXYBkg", "2-prong candidates (checked);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}});
      registry.add("MC/D0/hNormalisedDecLengthBkg", "2-prong candidates (checked);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}});
      registry.add("MC/D0/hNormalisedDecLengthXYBkg", "2-prong candidates (checked);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}});
      registry.add("MC/D0/hd0Prong0Bkg", "2-prong candidates (checked);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}});
      registry.add("MC/D0/hd0Prong1Bkg", "2-prong candidates (checked);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}});
      registry.add("MC/D0/hd0d0Bkg", "2-prong candidates (checked);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}});
      registry.add("MC/D0/hCTSBkg", "2-prong candidates (checked);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}});
      registry.add("MC/D0/hCtBkg", "2-prong candidates (checked);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}});
      registry.add("MC/D0/hCPABkg", "2-prong candidates (checked);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}});
      registry.add("MC/D0/hCPAxyBkg", "2-prong candidates (checked);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}});
      registry.add("MC/D0/hMassBkgD0", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}});
      registry.add("MC/D0/hMassReflBkgD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}});
      registry.add("MC/D0/hMassSigBkgD0bar", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}});
      registry.add("MC/D0/hMassSigD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}});
      registry.add("MC/D0/hMassReflBkgD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}});
      registry.add("MC/D0/hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/D0/hPtVsYGen", "2-prong candidates (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/D0/hYGenPrompt", "MC particles (matched, prompt);#it{y}^{gen.};entries", {HistType::kTH1F, {{300, -1.5, 1.5}}});
      registry.add("MC/D0/hPtVsYGenPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/D0/hYGenNonPrompt", "MC particles (matched, non-prompt);#it{y}^{gen.};entries", {HistType::kTH1F, {{300, -1.5, 1.5}}});
      registry.add("MC/D0/hPtVsYGenNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}});
      registry.add("MC/D0/hMassBkgD0bar", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigRecoTopol", "2-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigRecoCand", "2-prong candidates (RecoCand - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSig_RecoPID", "2-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtRecSig", "2-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/D0/hPtVsYRecSigPromptRecoHFFlag", "2-prong candidates (RecoHFFlag - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigPromptRecoTopol", "2-prong candidates (RecoTopol - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigPromptRecoCand", "2-prong candidates (RecoCand - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigPromptRecoPID", "2-prong candidates (RecoPID - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigPromptReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtRecSigPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/D0/hPtVsYRecSigNonPromptRecoHFFlag", "2-prong candidates (RecoHFFlag - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigNonPromptRecoTopol", "2-prong candidates (RecoTopol - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigNonPromptRecoCand", "2-prong candidates (RecoCand - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigNonPromptRecoPID", "2-prong candidates (RecoPID - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtVsYRecSigNonPromptReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}});
      registry.add("MC/D0/hPtRecSigNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});

      axesD0.push_back(thnAxisPtB);
      axesD0.push_back(thnAxisOrigin);
      registry.add("MC/D0/hMassVsPtVsFlatVsYVsD0Type", "Thn for D0 candidates", HistType::kTHnSparseD, axesD0);
      registry.add("MC/D0/hD0Gen", "Thn for generated D0 from charm and beauty", HistType::kTHnSparseD, axesD0Gen);
    }

    if (doprocessMCLc) {
      registry.add("MC/Lc/hPtGenSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hMassRecSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2});", {HistType::kTH1F, {{600, 1.98, 2.58}}});
      registry.add("MC/Lc/hMassVsPtRecSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins}}});
      registry.add("MC/Lc/hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hPtProng0RecSig", ";prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hPtProng1RecSig", ";prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hPtProng2RecSig", ";prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hd0Prong0RecSig", ";prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("MC/Lc/hd0Prong1RecSig", ";prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("MC/Lc/hd0Prong2RecSig", ";prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("MC/Lc/hd0VsPtProng0RecSig", ";prong 0 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("MC/Lc/hd0VsPtProng1RecSig", ";prong 1 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("MC/Lc/hd0VsPtProng2RecSig", ";prong 2 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("MC/Lc/hDecLengthRecSig", ";decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}});
      registry.add("MC/Lc/hDecLengthVsPtRecSig", ";decay length (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});
      registry.add("MC/Lc/hDecLengthxyRecSig", ";decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}});
      registry.add("MC/Lc/hDecLengthxyVsPtRecSig", ";decay length xy (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});
      registry.add("MC/Lc/hCtRecSig", ";proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}});
      registry.add("MC/Lc/hCtVsPtRecSig", ";proper lifetime (#Lambda_{c}) * #it{c} (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 0.2}, {vbins}}});
      registry.add("MC/Lc/hCPARecSig", ";cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("MC/Lc/hCPAVsPtRecSig", ";cosine of pointing angle;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});
      registry.add("MC/Lc/hCPAxyRecSig", ";cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("MC/Lc/hCPAxyVsPtRecSig", ";cosine of pointing angle xy;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});
      registry.add("MC/Lc/hEtaRecSig", ";#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("MC/Lc/hEtaVsPtRecSig", ";candidate #it{#eta};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
      registry.add("MC/Lc/hPhiRecSig", ";#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}});
      registry.add("MC/Lc/hPhiVsPtRecSig", ";candidate #it{#Phi};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});
      registry.add("MC/Lc/hImpParErrProng0VsPtRecSig", ";prong 0 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("MC/Lc/hImpParErrProng1VsPtRecSig", ";prong 1 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("MC/Lc/hImpParErrProng2VsPtRecSig", ";prong 2 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("MC/Lc/hDecLenErrVsPtRecSig", ";decay length error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 1.}, {vbins}}});
      registry.add("MC/Lc/hPtRecSigPrompt", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hPtProng0RecSigPrompt", ";prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hPtProng1RecSigPrompt", ";prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hPtProng2RecSigPrompt", ";prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hd0Prong0RecSigPrompt", ";prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("MC/Lc/hd0Prong1RecSigPrompt", ";prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("MC/Lc/hd0Prong2RecSigPrompt", ";prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("MC/Lc/hd0VsPtProng0RecSigPrompt", ";prong 0 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("MC/Lc/hd0VsPtProng1RecSigPrompt", ";prong 1 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("MC/Lc/hd0VsPtProng2RecSigPrompt", ";prong 2 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("MC/Lc/hDecLengthRecSigPrompt", ";decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}});
      registry.add("MC/Lc/hDecLengthVsPtRecSigPrompt", ";decay length (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});
      registry.add("MC/Lc/hDecLengthxyRecSigPrompt", ";decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}});
      registry.add("MC/Lc/hDecLengthxyVsPtRecSigPrompt", ";decay length xy (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});
      registry.add("MC/Lc/hCtRecSigPrompt", ";proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}});
      registry.add("MC/Lc/hCtVsPtRecSigPrompt", ";proper lifetime (#Lambda_{c}) * #it{c} (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 0.2}, {vbins}}});
      registry.add("MC/Lc/hCPARecSigPrompt", ";cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("MC/Lc/hCPAVsPtRecSigPrompt", ";cosine of pointing angle;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});
      registry.add("MC/Lc/hCPAxyRecSigPrompt", ";cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("MC/Lc/hCPAxyVsPtRecSigPrompt", ";cosine of pointing angle xy;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});
      registry.add("MC/Lc/hEtaRecSigPrompt", ";#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("MC/Lc/hEtaVsPtRecSigPrompt", ";candidate #it{#eta};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
      registry.add("MC/Lc/hPhiRecSigPrompt", ";#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}});
      registry.add("MC/Lc/hPhiVsPtRecSigPrompt", ";candidate #it{#Phi};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});
      registry.add("MC/Lc/hImpParErrProng0VsPtRecSigPrompt", ";prong 0 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("MC/Lc/hImpParErrProng1VsPtRecSigPrompt", ";prong 1 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("MC/Lc/hImpParErrProng2VsPtRecSigPrompt", ";prong 2 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("MC/Lc/hDecLenErrVsPtRecSigPrompt", ";decay length error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 1.}, {vbins}}});
      registry.add("MC/Lc/hPtRecSigNonPrompt", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hPtProng0RecSigNonPrompt", ";prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hPtProng1RecSigNonPrompt", ";prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hPtProng2RecSigNonPrompt", ";prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hd0Prong0RecSigNonPrompt", ";prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("MC/Lc/hd0Prong1RecSigNonPrompt", ";prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("MC/Lc/hd0Prong2RecSigNonPrompt", ";prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
      registry.add("MC/Lc/hd0VsPtProng0RecSigNonPrompt", ";prong 0 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("MC/Lc/hd0VsPtProng1RecSigNonPrompt", ";prong 1 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("MC/Lc/hd0VsPtProng2RecSigNonPrompt", ";prong 2 DCAxy to prim. vertex (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
      registry.add("MC/Lc/hDecLengthRecSigNonPrompt", ";decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}});
      registry.add("MC/Lc/hDecLengthVsPtRecSigNonPrompt", ";decay length (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});
      registry.add("MC/Lc/hDecLengthxyRecSigNonPrompt", ";decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}});
      registry.add("MC/Lc/hDecLengthxyVsPtRecSigNonPrompt", ";decay length xy (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});
      registry.add("MC/Lc/hCtRecSigNonPrompt", ";proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}});
      registry.add("MC/Lc/hCtVsPtRecSigNonPrompt", ";proper lifetime (#Lambda_{c}) * #it{c} (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 0.2}, {vbins}}});
      registry.add("MC/Lc/hCPARecSigNonPrompt", ";cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("MC/Lc/hCPAVsPtRecSigNonPrompt", ";cosine of pointing angle;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});
      registry.add("MC/Lc/hCPAxyRecSigNonPrompt", ";cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("MC/Lc/hCPAxyVsPtRecSigNonPrompt", ";cosine of pointing angle xy;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});
      registry.add("MC/Lc/hEtaRecSigNonPrompt", ";#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("MC/Lc/hEtaVsPtRecSigNonPrompt", ";candidate #it{#eta};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
      registry.add("MC/Lc/hPhiRecSigNonPrompt", ";#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}});
      registry.add("MC/Lc/hPhiVsPtRecSigNonPrompt", ";candidate #it{#Phi};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});
      registry.add("MC/Lc/hImpParErrProng0VsPtRecSigNonPrompt", ";prong 0 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("MC/Lc/hImpParErrProng1VsPtRecSigNonPrompt", ";prong 1 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("MC/Lc/hImpParErrProng2VsPtRecSigNonPrompt", ";prong 2 impact parameter error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
      registry.add("MC/Lc/hDecLenErrVsPtRecSigNonPrompt", ";decay length error (cm);#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 1.}, {vbins}}});
      registry.add("MC/Lc/hPtGen", ";#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hEtaGen", ";#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("MC/Lc/hYGen", ";#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("MC/Lc/hPhiGen", ";#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}});
      registry.add("MC/Lc/hEtaVsPtGen", ";#it{#eta};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
      registry.add("MC/Lc/hYVsPtGen", ";#it{y};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
      registry.add("MC/Lc/hPhiVsPtGen", ";#it{#Phi};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});
      registry.add("MC/Lc/hPtGenPrompt", ";#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hEtaGenPrompt", ";#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("MC/Lc/hYGenPrompt", ";#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("MC/Lc/hPhiGenPrompt", ";#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}});
      registry.add("MC/Lc/hEtaVsPtGenPrompt", ";#it{#eta};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
      registry.add("MC/Lc/hYVsPtGenPrompt", ";#it{y};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
      registry.add("MC/Lc/hPhiVsPtGenPrompt", ";#it{#Phi};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});
      registry.add("MC/Lc/hPtGenNonPrompt", ";#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("MC/Lc/hEtaGenNonPrompt", ";#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("MC/Lc/hYGenNonPrompt", ";#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("MC/Lc/hPhiGenNonPrompt", ";#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}});
      registry.add("MC/Lc/hEtaVsPtGenNonPrompt", ";#it{#eta};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
      registry.add("MC/Lc/hYVsPtGenNonPrompt", ";#it{y};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
      registry.add("MC/Lc/hPhiVsPtGenNonPrompt", ";#it{#Phi};#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});

      registry.add("MC/Lc/hnLcVars", "THn for Reconstructed Lambdac candidates for MC", HistType::kTHnSparseF, axesLcMc);
      registry.add("MC/Lc/hnLcVarsGen", "THn for Generated Lambdac", HistType::kTHnSparseF, axesGenLc);
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void processData(Collisions const& collisions,
                   D0Candidates const&,
                   LcCandidates const& selectedLcCandidates,
                   aod::TracksWExtra const&,
                   aod::BcFullInfos const& bcs,
                   aod::FV0As const&,
                   TracksWPid const& tracks)
  {
    runAnalysisData(collisions, selectedD0Candidates, selectedLcCandidates, bcs, tracks);
  }
  PROCESS_SWITCH(FlattenicityDLc, processData, "Process data", false);

  void processMCD0(D0CandidatesMc const&,
                   soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles2prong,
                   TracksSelQuality const& tracks,
                   CollisionsWithMcLabels const& collisions,
                   aod::McCollisions const& mcCollisions,
                   aod::FV0As const&,
                   aod::BcFullInfos const& bcs)
  {
    runAnalysisMCD0<aod::hf_cand::VertexerType::DCAFitter>(selectedD0CandidatesMc, mcParticles2prong, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(FlattenicityDLc, processMCD0, "Process MC D0 with DCAFitter", true);

  void processMCLc(soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticles3prong,
                   TracksSelQuality const& tracks,
                   CollisionsWithMcLabels const& collisions,
                   aod::McCollisions const& mcCollisions,
                   aod::FV0As const&,
                   LcCandidatesMc const& selectedLcCandidatesMc,
                   aod::BcFullInfos const& bcs)
  {
    runAnalysisMCLc<aod::hf_cand::VertexerType::DCAFitter>(selectedLcCandidatesMc, mcParticles3prong, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(FlattenicityDLc, processMCLc, "Process MC Lc with DCAFitter", true);

  template <typename CollType, typename CandTypeD0, typename CandTypeLc, typename BCsType>
  void runAnalysisData(CollType const& collisions,
                       CandTypeD0 const& candidatesD0,
                       CandTypeLc const& candidatesLc,
                       BCsType const& bcs,
                       TracksWPid const& tracks)
  {
    for (const auto& collision : collisions) {

      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, BCsType>(collision, centrality, ccdb, registry);
      hfEvSel.fillHistograms(collision, rejectionMask, centrality);
      if (rejectionMask != 0) {
        continue;
      }

      const float flat = fillFlat<true>(collision, 0);
      const float flat_calibrated = fillFlat<true>(collision, 1);

      const auto thisCollId = collision.globalIndex();

      // D0
      const auto& groupedD0Candidates = candidatesD0.sliceBy(candD0PerCollision, thisCollId);
      for (const auto& candidate : groupedD0Candidates) {
        if (yCandRecoMax >= 0. && std::abs(HfHelper::yD0(candidate)) > yCandRecoMax) {
          continue;
        }

        const float massD0 = HfHelper::invMassD0ToPiK(candidate);
        const float massD0bar = HfHelper::invMassD0barToKPi(candidate);
        const auto ptCandidate = candidate.pt();

        registry.fill(HIST("Data/D0/hPtCand"), ptCandidate);
        registry.fill(HIST("Data/D0/hPtProng0"), candidate.ptProng0());
        registry.fill(HIST("Data/D0/hPtProng1"), candidate.ptProng1());
        registry.fill(HIST("Data/D0/hDecLength"), candidate.decayLength(), ptCandidate);
        registry.fill(HIST("Data/D0/hDecLengthxy"), candidate.decayLengthXY(), ptCandidate);
        registry.fill(HIST("Data/D0/hDecLenErr"), candidate.errorDecayLength(), ptCandidate);
        registry.fill(HIST("Data/D0/hDecLenXYErr"), candidate.errorDecayLengthXY(), ptCandidate);
        registry.fill(HIST("Data/D0/hNormalisedDecLength"), candidate.decayLengthNormalised(), ptCandidate);
        registry.fill(HIST("Data/D0/hNormalisedDecLengthxy"), candidate.decayLengthXYNormalised(), ptCandidate);
        registry.fill(HIST("Data/D0/hd0Prong0"), candidate.impactParameter0(), ptCandidate);
        registry.fill(HIST("Data/D0/hd0Prong1"), candidate.impactParameter1(), ptCandidate);
        registry.fill(HIST("Data/D0/hd0ErrProng0"), candidate.errorImpactParameter0(), ptCandidate);
        registry.fill(HIST("Data/D0/hd0ErrProng1"), candidate.errorImpactParameter1(), ptCandidate);
        registry.fill(HIST("Data/D0/hd0d0"), candidate.impactParameterProduct(), ptCandidate);
        registry.fill(HIST("Data/D0/hCTS"), HfHelper::cosThetaStarD0(candidate), ptCandidate);
        registry.fill(HIST("Data/D0/hCt"), HfHelper::ctD0(candidate), ptCandidate);
        registry.fill(HIST("Data/D0/hCPA"), candidate.cpa(), ptCandidate);
        registry.fill(HIST("Data/D0/hEta"), candidate.eta(), ptCandidate);
        registry.fill(HIST("Data/D0/hSelectionStatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), ptCandidate);
        registry.fill(HIST("Data/D0/hDecLengthFinerBinning"), candidate.decayLength(), ptCandidate);
        registry.fill(HIST("Data/D0/hDecLengthxyFinerBinning"), candidate.decayLengthXY(), ptCandidate);
        registry.fill(HIST("Data/D0/hd0Prong0FinerBinning"), candidate.impactParameter0(), ptCandidate);
        registry.fill(HIST("Data/D0/hd0Prong1FinerBinning"), candidate.impactParameter1(), ptCandidate);
        registry.fill(HIST("Data/D0/hd0d0FinerBinning"), candidate.impactParameterProduct(), ptCandidate);
        registry.fill(HIST("Data/D0/hCTSFinerBinning"), HfHelper::cosThetaStarD0(candidate), ptCandidate);
        registry.fill(HIST("Data/D0/hCtFinerBinning"), HfHelper::ctD0(candidate), ptCandidate);
        registry.fill(HIST("Data/D0/hCPAFinerBinning"), candidate.cpa(), ptCandidate);
        registry.fill(HIST("Data/D0/hCPAXYFinerBinning"), candidate.cpaXY(), ptCandidate);

        if (candidate.isSelD0() >= selectionFlagD0) {
          registry.fill(HIST("Data/D0/hMass"), massD0, ptCandidate);
          registry.fill(HIST("Data/D0/hMassFinerBinning"), massD0, ptCandidate);
          registry.fill(HIST("Data/D0/hMassVsPhi"), massD0, ptCandidate, candidate.phi());

          registry.fill(HIST("Data/D0/hMassVsPtVsFlatVsYVsD0Type"), massD0, ptCandidate, flat, HfHelper::yD0(candidate), SigD0);
          registry.fill(HIST("Data/D0/hMassVsPtVsFlatVsYVsD0Type"), massD0, ptCandidate, flat, HfHelper::yD0(candidate), candidate.isSelD0bar() ? ReflectedD0 : PureSigD0);
        }
        if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          registry.fill(HIST("Data/D0/hMass"), massD0bar, ptCandidate);
          registry.fill(HIST("Data/D0/hMassFinerBinning"), massD0bar, ptCandidate);
          registry.fill(HIST("Data/D0/hMassVsPhi"), massD0bar, ptCandidate, candidate.phi());

          registry.fill(HIST("Data/D0/hMassVsPtVsFlatVsYVsD0Type"), massD0bar, ptCandidate, flat, HfHelper::yD0(candidate), SigD0bar);
          registry.fill(HIST("Data/D0/hMassVsPtVsFlatVsYVsD0Type"), massD0bar, ptCandidate, flat, HfHelper::yD0(candidate), candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar);
        }
      }

      // Lc
      const auto& groupedLcCandidates = candidatesLc.sliceBy(candLcPerCollision, thisCollId);

      for (const auto& candidate : groupedLcCandidates) {
        if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
          continue;
        }
        if (yCandRecoMax >= 0. && std::abs(HfHelper::yLc(candidate)) > yCandRecoMax) {
          continue;
        }
        const auto pt = candidate.pt();
        const auto ptProng0 = candidate.ptProng0();
        const auto ptProng1 = candidate.ptProng1();
        const auto ptProng2 = candidate.ptProng2();
        const auto decayLength = candidate.decayLength();
        const auto decayLengthXY = candidate.decayLengthXY();
        const auto chi2PCA = candidate.chi2PCA();
        const auto cpa = candidate.cpa();
        const auto cpaXY = candidate.cpaXY();

        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          registry.fill(HIST("Data/Lc/hMass"), HfHelper::invMassLcToPKPi(candidate));
          registry.fill(HIST("Data/Lc/hMassVsPt"), HfHelper::invMassLcToPKPi(candidate), pt);
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          registry.fill(HIST("Data/Lc/hMass"), HfHelper::invMassLcToPiKP(candidate));
          registry.fill(HIST("Data/Lc/hMassVsPt"), HfHelper::invMassLcToPiKP(candidate), pt);
        }
        registry.fill(HIST("Data/Lc/hPt"), pt);
        registry.fill(HIST("Data/Lc/hPtProng0"), ptProng0);
        registry.fill(HIST("Data/Lc/hPtProng1"), ptProng1);
        registry.fill(HIST("Data/Lc/hPtProng2"), ptProng2);
        registry.fill(HIST("Data/Lc/hd0Prong0"), candidate.impactParameter0());
        registry.fill(HIST("Data/Lc/hd0Prong1"), candidate.impactParameter1());
        registry.fill(HIST("Data/Lc/hd0Prong2"), candidate.impactParameter2());
        registry.fill(HIST("Data/Lc/hd0VsPtProng0"), candidate.impactParameter0(), pt);
        registry.fill(HIST("Data/Lc/hd0VsPtProng1"), candidate.impactParameter1(), pt);
        registry.fill(HIST("Data/Lc/hd0VsPtProng2"), candidate.impactParameter2(), pt);
        registry.fill(HIST("Data/Lc/hDecLength"), decayLength);
        registry.fill(HIST("Data/Lc/hDecLengthVsPt"), decayLength, pt);
        registry.fill(HIST("Data/Lc/hDecLengthxy"), decayLengthXY);
        registry.fill(HIST("Data/Lc/hDecLengthxyVsPt"), decayLengthXY, pt);
        registry.fill(HIST("Data/Lc/hCt"), HfHelper::ctLc(candidate));
        registry.fill(HIST("Data/Lc/hCtVsPt"), HfHelper::ctLc(candidate), pt);
        registry.fill(HIST("Data/Lc/hCPA"), cpa);
        registry.fill(HIST("Data/Lc/hCPAVsPt"), cpa, pt);
        registry.fill(HIST("Data/Lc/hCPAxy"), cpaXY);
        registry.fill(HIST("Data/Lc/hCPAxyVsPt"), cpaXY, pt);
        registry.fill(HIST("Data/Lc/hEta"), candidate.eta());
        registry.fill(HIST("Data/Lc/hEtaVsPt"), candidate.eta(), pt);
        registry.fill(HIST("Data/Lc/hPhi"), candidate.phi());
        registry.fill(HIST("Data/Lc/hPhiVsPt"), candidate.phi(), pt);
        registry.fill(HIST("Data/Lc/hSelectionStatus"), candidate.isSelLcToPKPi(), pt);
        registry.fill(HIST("Data/Lc/hSelectionStatus"), candidate.isSelLcToPiKP(), pt);
        registry.fill(HIST("Data/Lc/hImpParErrProng0VsPt"), candidate.errorImpactParameter0(), pt);
        registry.fill(HIST("Data/Lc/hImpParErrProng1VsPt"), candidate.errorImpactParameter1(), pt);
        registry.fill(HIST("Data/Lc/hImpParErrProng2VsPt"), candidate.errorImpactParameter2(), pt);
        registry.fill(HIST("Data/Lc/hDecLenErrVsPt"), candidate.errorDecayLength(), pt);

        auto fillTHnData = [&](bool isPKPi) {
          const auto massLc = isPKPi ? HfHelper::invMassLcToPKPi(candidate) : HfHelper::invMassLcToPiKP(candidate);
          std::vector<double> valuesToFill;
          valuesToFill.reserve(registry.get<THnSparse>(HIST("Data/Lc/hLcVsPtVsFlat"))->GetNdimensions());
          valuesToFill.insert(valuesToFill.end(), {massLc, pt, flat});

          registry.get<THnSparse>(HIST("Data/Lc/hLcVsPtVsFlat"))->Fill(valuesToFill.data());
        };

        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          fillTHnData(true);
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          fillTHnData(false);
        }
      }
    }
  }

  template <int ReconstructionType, typename CandTypeD0, typename CollType, typename BCsType>
  void runAnalysisMCD0(CandTypeD0 const& candidatesD0,
                       soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles2prong,
                       TracksSelQuality const&,
                       CollType const& collisions,
                       aod::McCollisions const&,
                       BCsType const&)
  {
    // MC rec.
    for (const auto& candidate : candidatesD0) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yD0(candidate)) > yCandRecoMax) {
        continue;
      }

      auto collision = candidate.template collision_as<CollType>();

      const float flat = fillFlat<false>(collision, 0);
      const float flat_calibrated = fillFlat<false>(collision, 1);

      float massD0{0.f}, massD0bar{0.f};
      massD0 = HfHelper::invMassD0ToPiK(candidate);
      massD0bar = HfHelper::invMassD0barToKPi(candidate);

      auto trackPos = candidate.template prong0_as<TracksSelQuality>();
      auto trackNeg = candidate.template prong1_as<TracksSelQuality>();
      if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
        auto indexMother = RecoDecay::getMother(mcParticles2prong, trackPos.template mcParticle_as<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>(), o2::constants::physics::Pdg::kD0, true);
        auto particleMother = mcParticles2prong.rawIteratorAt(indexMother);
        auto ptGen = particleMother.pt();
        auto yGen = RecoDecay::y(particleMother.pVector(), o2::constants::physics::MassD0);
        registry.fill(HIST("MC/D0/hPtGenSig"), ptGen);
        auto ptRec = candidate.pt();
        auto yRec = HfHelper::yD0(candidate);
        if (candidate.isRecoHfFlag() >= selectionFlagHf) {
          registry.fill(HIST("MC/D0/hPtVsYRecSigRecoHFFlag"), ptRec, yRec);
          registry.fill(HIST("MC/D0/hPtGenVsPtRecSig"), ptGen, ptRec);
          registry.fill(HIST("MC/D0/hYGenVsYRecSig"), yGen, yRec);
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("MC/D0/hMassVsPtGenVsPtRecSig"), massD0, ptGen, ptRec);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("MC/D0/hMassVsPtGenVsPtRecSig"), massD0bar, ptGen, ptRec);
          }
        }

        if (candidate.isRecoTopol() >= selectionTopol) {
          registry.fill(HIST("MC/D0/hPtVsYRecSigRecoTopol"), ptRec, yRec);
        }
        if (candidate.isRecoCand() >= selectionCand) {
          registry.fill(HIST("MC/D0/hPtVsYRecSigRecoCand"), ptRec, yRec);
        }
        if (candidate.isRecoPid() >= selectionPid) {
          registry.fill(HIST("MC/D0/hPtVsYRecSig_RecoPID"), ptRec, yRec);
        }
        if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0bar) {
          registry.fill(HIST("MC/D0/hPtVsYRecSigReco"), ptRec, yRec);
          registry.fill(HIST("MC/D0/hPtRecSig"), ptRec);
        }
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          if (candidate.isRecoHfFlag() >= selectionFlagHf) {
            registry.fill(HIST("MC/D0/hPtVsYRecSigPromptRecoHFFlag"), ptRec, yRec);
          }
          if (candidate.isRecoTopol() >= selectionTopol) {
            registry.fill(HIST("MC/D0/hPtVsYRecSigPromptRecoTopol"), ptRec, yRec);
          }
          if (candidate.isRecoCand() >= selectionCand) {
            registry.fill(HIST("MC/D0/hPtVsYRecSigPromptRecoCand"), ptRec, yRec);
          }
          if (candidate.isRecoPid() >= selectionPid) {
            registry.fill(HIST("MC/D0/hPtVsYRecSigPromptRecoPID"), ptRec, yRec);
          }
          if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("MC/D0/hPtVsYRecSigPromptReco"), ptRec, yRec);
            registry.fill(HIST("MC/D0/hPtRecSigPrompt"), ptRec);
          }
        } else {
          if (candidate.isRecoHfFlag() >= selectionFlagHf) {
            registry.fill(HIST("MC/D0/hPtVsYRecSigNonPromptRecoHFFlag"), ptRec, yRec);
          }
          if (candidate.isRecoTopol() >= selectionTopol) {
            registry.fill(HIST("MC/D0/hPtVsYRecSigNonPromptRecoTopol"), ptRec, yRec);
          }
          if (candidate.isRecoCand() >= selectionCand) {
            registry.fill(HIST("MC/D0/hPtVsYRecSigNonPromptRecoCand"), ptRec, yRec);
          }
          if (candidate.isRecoPid() >= selectionPid) {
            registry.fill(HIST("MC/D0/hPtVsYRecSigNonPromptRecoPID"), ptRec, yRec);
          }
          if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("MC/D0/hPtVsYRecSigNonPromptReco"), ptRec, yRec);
            registry.fill(HIST("MC/D0/hPtRecSigNonPrompt"), ptRec);
          }
        }
        registry.fill(HIST("MC/D0/hCPARecSig"), candidate.cpa());
        registry.fill(HIST("MC/D0/hEtaRecSig"), candidate.eta());
      } else {
        registry.fill(HIST("MC/D0/hPtRecBg"), candidate.pt());
        registry.fill(HIST("MC/D0/hCPARecBg"), candidate.cpa());
        registry.fill(HIST("MC/D0/hEtaRecBg"), candidate.eta());
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
        registry.fill(HIST("MC/D0/hMassSigBkgD0"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          registry.fill(HIST("MC/D0/hPtProng0Sig"), ptProng0, rapidityCandidate);
          registry.fill(HIST("MC/D0/hPtProng1Sig"), ptProng1, rapidityCandidate);
          registry.fill(HIST("MC/D0/hDecLengthSig"), declengthCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hDecLengthXYSig"), declengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hNormalisedDecLengthSig"), normaliseddeclengthCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hNormalisedDecLengthXYSig"), normaliseddeclengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hd0Prong0Sig"), d0Prong0, rapidityCandidate);
          registry.fill(HIST("MC/D0/hd0Prong1Sig"), d0Prong1, rapidityCandidate);
          registry.fill(HIST("MC/D0/hd0d0Sig"), d0d0Candidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hCTSSig"), ctsCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hCtSig"), ctCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hCPASig"), cpaCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hCPAxySig"), cpaxyCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hd0Prong0VsPtSig"), d0Prong0, ptCandidate);
          registry.fill(HIST("MC/D0/hd0Prong1VsPtSig"), d0Prong1, ptCandidate);
          registry.fill(HIST("MC/D0/hd0d0VsPtSig"), d0d0Candidate, ptCandidate);
          registry.fill(HIST("MC/D0/hCTSVsPtSig"), ctsCandidate, ptCandidate);
          registry.fill(HIST("MC/D0/hCPAVsPtSig"), cpaCandidate, ptCandidate);
          registry.fill(HIST("MC/D0/hCPAXYVsPtSig"), cpaxyCandidate, ptCandidate);
          registry.fill(HIST("MC/D0/hNormalisedDecLengthxyVsPtSig"), normaliseddeclengthxyCandidate, ptCandidate);
          registry.fill(HIST("MC/D0/hDecLengthVsPtSig"), declengthCandidate, ptCandidate);
          registry.fill(HIST("MC/D0/hDecLengthxyVsPtSig"), declengthxyCandidate, ptCandidate);
          registry.fill(HIST("MC/D0/hMassSigD0"), massD0, ptCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hMassVsPtVsFlatVsYVsD0Type"), massD0, ptCandidate, flat, rapidityCandidate, SigD0, candidate.ptBhadMotherPart(), candidate.originMcRec());
        } else {
          registry.fill(HIST("MC/D0/hPtProng0Bkg"), ptProng0, rapidityCandidate);
          registry.fill(HIST("MC/D0/hPtProng1Bkg"), ptProng1, rapidityCandidate);
          registry.fill(HIST("MC/D0/hDecLengthBkg"), declengthCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hDecLengthXYBkg"), declengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hNormalisedDecLengthBkg"), normaliseddeclengthCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hNormalisedDecLengthXYBkg"), normaliseddeclengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hd0Prong0Bkg"), d0Prong0, rapidityCandidate);
          registry.fill(HIST("MC/D0/hd0Prong1Bkg"), d0Prong1, rapidityCandidate);
          registry.fill(HIST("MC/D0/hd0d0Bkg"), d0d0Candidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hCTSBkg"), ctsCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hCtBkg"), ctCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hCPABkg"), cpaCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hCPAxyBkg"), cpaxyCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hMassBkgD0"), massD0, ptCandidate, rapidityCandidate);
          if (candidate.flagMcMatchRec() == -o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
            registry.fill(HIST("MC/D0/hMassReflBkgD0"), massD0, ptCandidate, rapidityCandidate);
            registry.fill(HIST("MC/D0/hMassVsPtVsFlatVsYVsD0Type"), massD0, ptCandidate, flat, rapidityCandidate, ReflectedD0, candidate.ptBhadMotherPart(), candidate.originMcRec());
          }
        }
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("MC/D0/hMassSigBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == -o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          registry.fill(HIST("MC/D0/hMassSigD0bar"), massD0bar, ptCandidate, rapidityCandidate);
          registry.fill(HIST("MC/D0/hMassVsPtVsFlatVsYVsD0Type"), massD0bar, ptCandidate, flat, rapidityCandidate, SigD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec());
        } else {
          registry.fill(HIST("MC/D0/hMassBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
          if (candidate.flagMcMatchRec() == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
            registry.fill(HIST("MC/D0/hMassReflBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
            registry.fill(HIST("MC/D0/hMassVsPtVsFlatVsYVsD0Type"), massD0bar, ptCandidate, flat, rapidityCandidate, ReflectedD0bar, candidate.ptBhadMotherPart(), candidate.originMcRec());
          }
        }
      }
    }

    // MC gen.
    for (const auto& particle : mcParticles2prong) {
      if (std::abs(particle.flagMcMatchGen()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
        if (yCandGenMax >= 0. && std::abs(RecoDecay::y(particle.pVector(), o2::constants::physics::MassD0)) > yCandGenMax) {
          continue;
        }

        float flat{-1.f};
        const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, particle.mcCollision().globalIndex());
        for (const auto& recCol : recoCollsPerMcColl) {
          flat = fillFlat<false>(recCol, 0);
        }

        float ptGenB = -1;
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassD0);
        registry.fill(HIST("MC/D0/hPtGen"), ptGen);
        registry.fill(HIST("MC/D0/hPtVsYGen"), ptGen, yGen);

        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("MC/D0/hPtGenPrompt"), ptGen);
          registry.fill(HIST("MC/D0/hYGenPrompt"), yGen);
          registry.fill(HIST("MC/D0/hPtVsYGenPrompt"), ptGen, yGen);
          registry.fill(HIST("MC/D0/hD0Gen"), ptGen, flat, ptGenB, yGen, 1);
        } else {
          ptGenB = mcParticles2prong.rawIteratorAt(particle.idxBhadMotherPart()).pt();
          registry.fill(HIST("MC/D0/hPtGenNonPrompt"), ptGen);
          registry.fill(HIST("MC/D0/hYGenNonPrompt"), yGen);
          registry.fill(HIST("MC/D0/hPtVsYGenNonPrompt"), ptGen, yGen);
          registry.fill(HIST("MC/D0/hD0Gen"), ptGen, flat, ptGenB, yGen, 2);
        }
        registry.fill(HIST("MC/D0/hEtaGen"), particle.eta());
      }
    }
  }

  template <int ReconstructionType, typename CandTypeLc, typename CollType, typename BCsType>
  void runAnalysisMCLc(CandTypeLc const& candidatesLc,
                       soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticles3prong,
                       TracksSelQuality const&,
                       CollType const& collisions,
                       aod::McCollisions const&,
                       BCsType const&)
  {
    for (const auto& collision : collisions) {
      // MC Rec.
      const auto thisCollId = collision.globalIndex();
      const auto& groupedLcCandidates = candidatesLc.sliceBy(candLcPerCollision, thisCollId);

      const float flat = fillFlat<true>(collision, 0);
      const float flat_calibrated = fillFlat<true>(collision, 1);

      for (const auto& candidate : groupedLcCandidates) {
        if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
          continue;
        }
        if (yCandRecoMax >= 0. && std::abs(HfHelper::yLc(candidate)) > yCandRecoMax) {
          continue;
        }

        if (std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
          const auto& mcParticleProng0 = candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>();
          const auto pdgCodeProng0 = std::abs(mcParticleProng0.pdgCode());
          const auto indexMother = RecoDecay::getMother(mcParticles3prong, mcParticleProng0, o2::constants::physics::Pdg::kLambdaCPlus, true);
          const auto particleMother = mcParticles3prong.rawIteratorAt(indexMother);
          registry.fill(HIST("MC/Lc/hPtGenSig"), particleMother.pt());

          const auto pt = candidate.pt();
          const auto ptProng0 = candidate.ptProng0();
          const auto ptProng1 = candidate.ptProng1();
          const auto ptProng2 = candidate.ptProng2();
          const auto decayLength = candidate.decayLength();
          const auto chi2PCA = candidate.chi2PCA();
          const auto cpa = candidate.cpa();
          const auto originType = candidate.originMcRec();
          const auto ptRecB = candidate.ptBhadMotherPart();

          if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
            registry.fill(HIST("MC/Lc/hMassRecSig"), HfHelper::invMassLcToPKPi(candidate));
            registry.fill(HIST("MC/Lc/hMassVsPtRecSig"), HfHelper::invMassLcToPKPi(candidate), candidate.pt());
          }
          if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
            registry.fill(HIST("MC/Lc/hMassRecSig"), HfHelper::invMassLcToPiKP(candidate));
            registry.fill(HIST("MC/Lc/hMassVsPtRecSig"), HfHelper::invMassLcToPiKP(candidate), candidate.pt());
          }
          registry.fill(HIST("MC/Lc/hPtRecSig"), candidate.pt());
          registry.fill(HIST("MC/Lc/hPtProng0RecSig"), candidate.ptProng0());
          registry.fill(HIST("MC/Lc/hPtProng1RecSig"), candidate.ptProng1());
          registry.fill(HIST("MC/Lc/hPtProng2RecSig"), candidate.ptProng2());
          registry.fill(HIST("MC/Lc/hd0Prong0RecSig"), candidate.impactParameter0());
          registry.fill(HIST("MC/Lc/hd0Prong1RecSig"), candidate.impactParameter1());
          registry.fill(HIST("MC/Lc/hd0Prong2RecSig"), candidate.impactParameter2());
          registry.fill(HIST("MC/Lc/hd0VsPtProng0RecSig"), candidate.impactParameter0(), candidate.pt());
          registry.fill(HIST("MC/Lc/hd0VsPtProng1RecSig"), candidate.impactParameter1(), candidate.pt());
          registry.fill(HIST("MC/Lc/hd0VsPtProng2RecSig"), candidate.impactParameter2(), candidate.pt());
          registry.fill(HIST("MC/Lc/hDecLengthRecSig"), candidate.decayLength());
          registry.fill(HIST("MC/Lc/hDecLengthVsPtRecSig"), candidate.decayLength(), candidate.pt());
          registry.fill(HIST("MC/Lc/hDecLengthxyRecSig"), candidate.decayLengthXY());
          registry.fill(HIST("MC/Lc/hDecLengthxyVsPtRecSig"), candidate.decayLengthXY(), candidate.pt());
          registry.fill(HIST("MC/Lc/hCtRecSig"), HfHelper::ctLc(candidate));
          registry.fill(HIST("MC/Lc/hCtVsPtRecSig"), HfHelper::ctLc(candidate), candidate.pt());
          registry.fill(HIST("MC/Lc/hCPARecSig"), candidate.cpa());
          registry.fill(HIST("MC/Lc/hCPAVsPtRecSig"), candidate.cpa(), candidate.pt());
          registry.fill(HIST("MC/Lc/hCPAxyRecSig"), candidate.cpaXY());
          registry.fill(HIST("MC/Lc/hCPAxyVsPtRecSig"), candidate.cpaXY(), candidate.pt());
          registry.fill(HIST("MC/Lc/hEtaRecSig"), candidate.eta());
          registry.fill(HIST("MC/Lc/hEtaVsPtRecSig"), candidate.eta(), candidate.pt());
          registry.fill(HIST("MC/Lc/hPhiRecSig"), candidate.phi());
          registry.fill(HIST("MC/Lc/hPhiVsPtRecSig"), candidate.phi(), candidate.pt());
          registry.fill(HIST("MC/Lc/hImpParErrProng0VsPtRecSig"), candidate.errorImpactParameter0(), candidate.pt());
          registry.fill(HIST("MC/Lc/hImpParErrProng1VsPtRecSig"), candidate.errorImpactParameter1(), candidate.pt());
          registry.fill(HIST("MC/Lc/hImpParErrProng2VsPtRecSig"), candidate.errorImpactParameter2(), candidate.pt());
          registry.fill(HIST("MC/Lc/hDecLenErrVsPtRecSig"), candidate.errorDecayLength(), candidate.pt());

          if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
            registry.fill(HIST("MC/Lc/hPtRecSigPrompt"), candidate.pt());
            registry.fill(HIST("MC/Lc/hPtProng0RecSigPrompt"), candidate.ptProng0());
            registry.fill(HIST("MC/Lc/hPtProng1RecSigPrompt"), candidate.ptProng1());
            registry.fill(HIST("MC/Lc/hPtProng2RecSigPrompt"), candidate.ptProng2());
            registry.fill(HIST("MC/Lc/hd0Prong0RecSigPrompt"), candidate.impactParameter0());
            registry.fill(HIST("MC/Lc/hd0Prong1RecSigPrompt"), candidate.impactParameter1());
            registry.fill(HIST("MC/Lc/hd0Prong2RecSigPrompt"), candidate.impactParameter2());
            registry.fill(HIST("MC/Lc/hd0VsPtProng0RecSigPrompt"), candidate.impactParameter0(), candidate.pt());
            registry.fill(HIST("MC/Lc/hd0VsPtProng1RecSigPrompt"), candidate.impactParameter1(), candidate.pt());
            registry.fill(HIST("MC/Lc/hd0VsPtProng2RecSigPrompt"), candidate.impactParameter2(), candidate.pt());
            registry.fill(HIST("MC/Lc/hDecLengthRecSigPrompt"), candidate.decayLength());
            registry.fill(HIST("MC/Lc/hDecLengthVsPtRecSigPrompt"), candidate.decayLength(), candidate.pt());
            registry.fill(HIST("MC/Lc/hDecLengthxyRecSigPrompt"), candidate.decayLengthXY());
            registry.fill(HIST("MC/Lc/hDecLengthxyVsPtRecSigPrompt"), candidate.decayLengthXY(), candidate.pt());
            registry.fill(HIST("MC/Lc/hCtRecSigPrompt"), HfHelper::ctLc(candidate));
            registry.fill(HIST("MC/Lc/hCtVsPtRecSigPrompt"), HfHelper::ctLc(candidate), candidate.pt());
            registry.fill(HIST("MC/Lc/hCPARecSigPrompt"), candidate.cpa());
            registry.fill(HIST("MC/Lc/hCPAVsPtRecSigPrompt"), candidate.cpa(), candidate.pt());
            registry.fill(HIST("MC/Lc/hCPAxyRecSigPrompt"), candidate.cpaXY());
            registry.fill(HIST("MC/Lc/hCPAxyVsPtRecSigPrompt"), candidate.cpaXY(), candidate.pt());
            registry.fill(HIST("MC/Lc/hEtaRecSigPrompt"), candidate.eta());
            registry.fill(HIST("MC/Lc/hEtaVsPtRecSigPrompt"), candidate.eta(), candidate.pt());
            registry.fill(HIST("MC/Lc/hPhiRecSigPrompt"), candidate.phi());
            registry.fill(HIST("MC/Lc/hPhiVsPtRecSigPrompt"), candidate.phi(), candidate.pt());
            registry.fill(HIST("MC/Lc/hImpParErrProng0VsPtRecSigPrompt"), candidate.errorImpactParameter0(), candidate.pt());
            registry.fill(HIST("MC/Lc/hImpParErrProng1VsPtRecSigPrompt"), candidate.errorImpactParameter1(), candidate.pt());
            registry.fill(HIST("MC/Lc/hImpParErrProng2VsPtRecSigPrompt"), candidate.errorImpactParameter2(), candidate.pt());
            registry.fill(HIST("MC/Lc/hDecLenErrVsPtRecSigPrompt"), candidate.errorDecayLength(), candidate.pt());
          } else if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
            registry.fill(HIST("MC/Lc/hPtRecSigNonPrompt"), candidate.pt());
            registry.fill(HIST("MC/Lc/hPtProng0RecSigNonPrompt"), candidate.ptProng0());
            registry.fill(HIST("MC/Lc/hPtProng1RecSigNonPrompt"), candidate.ptProng1());
            registry.fill(HIST("MC/Lc/hPtProng2RecSigNonPrompt"), candidate.ptProng2());
            registry.fill(HIST("MC/Lc/hd0Prong0RecSigNonPrompt"), candidate.impactParameter0());
            registry.fill(HIST("MC/Lc/hd0Prong1RecSigNonPrompt"), candidate.impactParameter1());
            registry.fill(HIST("MC/Lc/hd0Prong2RecSigNonPrompt"), candidate.impactParameter2());
            registry.fill(HIST("MC/Lc/hd0VsPtProng0RecSigNonPrompt"), candidate.impactParameter0(), candidate.pt());
            registry.fill(HIST("MC/Lc/hd0VsPtProng1RecSigNonPrompt"), candidate.impactParameter1(), candidate.pt());
            registry.fill(HIST("MC/Lc/hd0VsPtProng2RecSigNonPrompt"), candidate.impactParameter2(), candidate.pt());
            registry.fill(HIST("MC/Lc/hDecLengthRecSigNonPrompt"), candidate.decayLength());
            registry.fill(HIST("MC/Lc/hDecLengthVsPtRecSigNonPrompt"), candidate.decayLength(), candidate.pt());
            registry.fill(HIST("MC/Lc/hDecLengthxyRecSigNonPrompt"), candidate.decayLengthXY());
            registry.fill(HIST("MC/Lc/hDecLengthxyVsPtRecSigNonPrompt"), candidate.decayLengthXY(), candidate.pt());
            registry.fill(HIST("MC/Lc/hCtRecSigNonPrompt"), HfHelper::ctLc(candidate));
            registry.fill(HIST("MC/Lc/hCtVsPtRecSigNonPrompt"), HfHelper::ctLc(candidate), candidate.pt());
            registry.fill(HIST("MC/Lc/hCPARecSigNonPrompt"), candidate.cpa());
            registry.fill(HIST("MC/Lc/hCPAVsPtRecSigNonPrompt"), candidate.cpa(), candidate.pt());
            registry.fill(HIST("MC/Lc/hCPAxyRecSigNonPrompt"), candidate.cpaXY());
            registry.fill(HIST("MC/Lc/hCPAxyVsPtRecSigNonPrompt"), candidate.cpaXY(), candidate.pt());
            registry.fill(HIST("MC/Lc/hEtaRecSigNonPrompt"), candidate.eta());
            registry.fill(HIST("MC/Lc/hEtaVsPtRecSigNonPrompt"), candidate.eta(), candidate.pt());
            registry.fill(HIST("MC/Lc/hPhiRecSigNonPrompt"), candidate.phi());
            registry.fill(HIST("MC/Lc/hPhiVsPtRecSigNonPrompt"), candidate.phi(), candidate.pt());
            registry.fill(HIST("MC/Lc/hImpParErrProng0VsPtRecSigNonPrompt"), candidate.errorImpactParameter0(), candidate.pt());
            registry.fill(HIST("MC/Lc/hImpParErrProng1VsPtRecSigNonPrompt"), candidate.errorImpactParameter1(), candidate.pt());
            registry.fill(HIST("MC/Lc/hImpParErrProng2VsPtRecSigNonPrompt"), candidate.errorImpactParameter2(), candidate.pt());
            registry.fill(HIST("MC/Lc/hDecLenErrVsPtRecSigNonPrompt"), candidate.errorDecayLength(), candidate.pt());
          }

          if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
            const auto massLc = HfHelper::invMassLcToPKPi(candidate);
            std::vector<double> valuesToFill;
            valuesToFill.reserve(registry.get<THnSparse>(HIST("MC/Lc/hnLcVars"))->GetNdimensions());
            valuesToFill.insert(valuesToFill.end(), {massLc, pt, flat, ptRecB, static_cast<double>(originType)});
            registry.get<THnSparse>(HIST("MC/Lc/hnLcVars"))->Fill(valuesToFill.data());
          }
          if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
            const auto massLc = HfHelper::invMassLcToPiKP(candidate);
            std::vector<double> valuesToFill;
            valuesToFill.reserve(registry.get<THnSparse>(HIST("MC/Lc/hnLcVars"))->GetNdimensions());
            valuesToFill.insert(valuesToFill.end(), {massLc, pt, flat, ptRecB, static_cast<double>(originType)});
            registry.get<THnSparse>(HIST("MC/Lc/hnLcVars"))->Fill(valuesToFill.data());
          }
        }
      }
    }

    // MC gen.
    for (const auto& particle : mcParticles3prong) {
      if (std::abs(particle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus);
        if (yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) {
          continue;
        }

        float flat{-1.f};
        const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollisionLc, particle.mcCollision().globalIndex());
        for (const auto& recCol : recoCollsPerMcColl) {
          flat = fillFlat<false>(recCol, 0);
        }

        const auto ptGen = particle.pt();
        const auto originType = particle.originMcGen();
        float ptGenB = -1.;

        registry.fill(HIST("MC/Lc/hPtGen"), particle.pt());
        registry.fill(HIST("MC/Lc/hEtaGen"), particle.eta());
        registry.fill(HIST("MC/Lc/hYGen"), RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus));
        registry.fill(HIST("MC/Lc/hPhiGen"), particle.phi());
        registry.fill(HIST("MC/Lc/hEtaVsPtGen"), particle.eta(), particle.pt());
        registry.fill(HIST("MC/Lc/hYVsPtGen"), RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus), particle.pt());
        registry.fill(HIST("MC/Lc/hPhiVsPtGen"), particle.phi(), particle.pt());

        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          ptGenB = -1.;
          std::vector<double> valuesToFill{ptGen, flat, ptGenB, static_cast<double>(originType)};
          registry.get<THnSparse>(HIST("MC/Lc/hnLcVarsGen"))->Fill(valuesToFill.data());

          registry.fill(HIST("MC/Lc/hPtGenPrompt"), particle.pt());
          registry.fill(HIST("MC/Lc/hEtaGenPrompt"), particle.eta());
          registry.fill(HIST("MC/Lc/hYGenPrompt"), RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus));
          registry.fill(HIST("MC/Lc/hPhiGenPrompt"), particle.phi());
          registry.fill(HIST("MC/Lc/hEtaVsPtGenPrompt"), particle.eta(), particle.pt());
          registry.fill(HIST("MC/Lc/hYVsPtGenPrompt"), RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus), particle.pt());
          registry.fill(HIST("MC/Lc/hPhiVsPtGenPrompt"), particle.phi(), particle.pt());
        } else if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          ptGenB = mcParticles3prong.rawIteratorAt(particle.idxBhadMotherPart()).pt();
          std::vector<double> valuesToFill{ptGen, flat, ptGenB, static_cast<double>(originType)};
          registry.get<THnSparse>(HIST("MC/Lc/hnLcVarsGen"))->Fill(valuesToFill.data());

          registry.fill(HIST("MC/Lc/hPtGenNonPrompt"), particle.pt());
          registry.fill(HIST("MC/Lc/hEtaGenNonPrompt"), particle.eta());
          registry.fill(HIST("MC/Lc/hYGenNonPrompt"), RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus));
          registry.fill(HIST("MC/Lc/hPhiGenNonPrompt"), particle.phi());
          registry.fill(HIST("MC/Lc/hEtaVsPtGenNonPrompt"), particle.eta(), particle.pt());
          registry.fill(HIST("MC/Lc/hYVsPtGenNonPrompt"), RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus), particle.pt());
          registry.fill(HIST("MC/Lc/hPhiVsPtGenNonPrompt"), particle.phi(), particle.pt());
        }
      }
    }
  }

  template <bool fillHist = true, typename CollType>
  float fillFlat(CollType const& collision, bool const& ifCalib)
  {
    rhoLatticeFV0.fill(0);
    fv0AmplitudeWoCalib.fill(0);
    if (collision.has_foundFV0()) {
      auto fv0 = collision.foundFV0();
      std::bitset<8> fV0Triggers = fv0.triggerMask();
      bool isOkFV0OrA = fV0Triggers[o2::fit::Triggers::bitA];
      if (isOkFV0OrA) {
        for (std::size_t ich = 0; ich < fv0.channel().size(); ich++) {
          float amplCh = fv0.amplitude()[ich];
          int chv0 = fv0.channel()[ich];
          int chv0phi = getFV0IndexPhi(chv0);
          if (amplCh > 0.0) {
            if (chv0phi > 0.0) {
              fv0AmplitudeWoCalib[chv0phi] = amplCh;
              if (ifCalib) {
                amplCh *= calib[chv0phi];
              }
              if (chv0 < CinnerFV0) {
                rhoLatticeFV0[chv0phi] += amplCh;
              } else {
                rhoLatticeFV0[chv0phi] += amplCh / 2.;
              }
            }
          }
        }
        float flattenicityFV0 = calcFlatenicity(rhoLatticeFV0);
        if constexpr (fillHist) {
          if (ifCalib)
            registry.fill(HIST("Flattenicity_calibrated"), 1 - flattenicityFV0);
          else
            registry.fill(HIST("Flattenicity"), 1 - flattenicityFV0);
        }
        return 1. - flattenicityFV0;
      } else {
        return 9999;
      }
    } else {
      return 9999;
    }
  }

  template <typename T, std::size_t S>
  float calcFlatenicity(std::array<T, S> const& signals)
  {
    static_assert(S != 0);

    int entries = signals.size();
    float flat{-1};
    float mRho{0};
    float mRho_debug{0};
    for (int iCell = 0; iCell < entries; ++iCell) {
      if (signals[iCell] > 0.0) {
        mRho += 1.0 * signals[iCell];
        mRho_debug += 1.0 * signals[iCell];
      }
    }
    mRho /= (1.0 * entries);
    if (mRho <= 0) {
      return -1;
    }
    float sRhoTmp{0};
    float sRho{0};
    for (int iCell = 0; iCell < entries; ++iCell) {
      if (signals[iCell] > 0.0) {
        sRhoTmp += std::pow(1.0 * signals[iCell] - mRho, 2);
      }
    }
    sRhoTmp /= (1.0 * entries * entries);
    sRho = std::sqrt(sRhoTmp);
    if (mRho > 0.0) {
      flat = sRho / mRho;
    } else {
      flat = -1;
    }
    return flat;
  }

  int getFV0IndexPhi(int i_ch)
  {
    int iRing = -1;

    if (i_ch >= 0 && i_ch < 8) {
      if (i_ch < 4) {
        iRing = i_ch;
      } else {
        if (i_ch == 7) {
          iRing = 4;
        } else if (i_ch == 6) {
          iRing = 5;
        } else if (i_ch == 5) {
          iRing = 6;
        } else if (i_ch == 4) {
          iRing = 7;
        }
      }
    } else if (i_ch >= 8 && i_ch < 16) {
      if (i_ch < 12) {
        iRing = i_ch;
      } else {
        if (i_ch == 15) {
          iRing = 12;
        } else if (i_ch == 14) {
          iRing = 13;
        } else if (i_ch == 13) {
          iRing = 14;
        } else if (i_ch == 12) {
          iRing = 15;
        }
      }
    } else if (i_ch >= 16 && i_ch < 24) {
      if (i_ch < 20) {
        iRing = i_ch;
      } else {
        if (i_ch == 23) {
          iRing = 20;
        } else if (i_ch == 22) {
          iRing = 21;
        } else if (i_ch == 21) {
          iRing = 22;
        } else if (i_ch == 20) {
          iRing = 23;
        }
      }
    } else if (i_ch >= 24 && i_ch < 32) {
      if (i_ch < 28) {
        iRing = i_ch;
      } else {
        if (i_ch == 31) {
          iRing = 28;
        } else if (i_ch == 30) {
          iRing = 29;
        } else if (i_ch == 29) {
          iRing = 30;
        } else if (i_ch == 28) {
          iRing = 31;
        }
      }
    } else if (i_ch == 32) {
      iRing = 32;
    } else if (i_ch == 40) {
      iRing = 33;
    } else if (i_ch == 33) {
      iRing = 34;
    } else if (i_ch == 41) {
      iRing = 35;
    } else if (i_ch == 34) {
      iRing = 36;
    } else if (i_ch == 42) {
      iRing = 37;
    } else if (i_ch == 35) {
      iRing = 38;
    } else if (i_ch == 43) {
      iRing = 39;
    } else if (i_ch == 47) {
      iRing = 40;
    } else if (i_ch == 39) {
      iRing = 41;
    } else if (i_ch == 46) {
      iRing = 42;
    } else if (i_ch == 38) {
      iRing = 43;
    } else if (i_ch == 45) {
      iRing = 44;
    } else if (i_ch == 37) {
      iRing = 45;
    } else if (i_ch == 44) {
      iRing = 46;
    } else if (i_ch == 36) {
      iRing = 47;
    }
    return iRing;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FlattenicityDLc>(cfgc)};
}
