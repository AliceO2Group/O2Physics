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

/// \file correlatorLcHadrons.cxx
/// \brief Lc-Hadrons correlator task - data-like, Mc-Reco and Mc-Gen analyses
///
/// \author Marianna Mazzilli <marianna.mazzilli@cern.ch>
/// \author Zhen Zhang <zhenz@cern.ch>
/// \author Ravindra Singh <ravindra.singh@cern.ch>

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFC/Utils/utilsCorrelations.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::hf_correlations;
///
/// Returns deltaPhi values in range [-pi/2., 3.*pi/2.], typically used for correlation studies
///
double getDeltaPhi(double phiLc, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiLc, -PIHalf);
}

// definition of ME variables
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
using BinningTypeMcGen = ColumnBinningPolicy<aod::mccollision::PosZ, o2::aod::mult::MultMCFT0A>;

// Code to select collisions with at least one Lambda_c
struct HfCorrelatorLcHadronsSelection {
  Produces<aod::LcSelection> lcSel;

  Configurable<bool> useSel8{"useSel8", true, "Flag for applying sel8 for collision selection"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<bool> doSelLcCollision{"doSelLcCollision", true, "Select collisions with at least one Lc"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  HfHelper hfHelper;
  SliceCache cache;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CandidatesLcData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using CandidatesLcMcRec = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>;
  using CandidatesLcMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  // filter on selection of Lc and decay channel Lc->PKPi
  Filter lcFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);

  /// Code to select collisions with at least one Lc - for real data and data-like analysis
  void processLcSelectionData(SelCollisions::iterator const& collision,
                              CandidatesLcData const& candidates)
  {
    bool isSelColl = true;
    bool isLcFound = true;
    bool isSel8 = true;
    bool isNosameBunchPileUp = true;
    if (doSelLcCollision) {
      for (const auto& candidate : candidates) {
        if (std::abs(hfHelper.yLc(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          isLcFound = false;
          continue;
        }
        isLcFound = true;
        break;
      }
    }
    if (useSel8) {
      isSel8 = collision.sel8();
    }
    if (selNoSameBunchPileUpColl) {
      isNosameBunchPileUp = static_cast<bool>(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup));
    }
    isSelColl = isLcFound && isSel8 && isNosameBunchPileUp;
    lcSel(isSelColl);
  }
  PROCESS_SWITCH(HfCorrelatorLcHadronsSelection, processLcSelectionData, "Process Lc Collision Selection Data", true);

  void processLcSelectionMcRec(SelCollisions::iterator const& collision,
                               CandidatesLcMcRec const& candidates)
  {
    bool isSelColl = true;
    bool isLcFound = true;
    bool isSel8 = true;
    bool isNosameBunchPileUp = true;
    if (doSelLcCollision) {
      for (const auto& candidate : candidates) {
        if (std::abs(hfHelper.yLc(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          isLcFound = false;
          continue;
        }
        isLcFound = true;
        break;
      }
    }
    if (useSel8) {
      isSel8 = collision.sel8();
    }
    if (selNoSameBunchPileUpColl) {
      isNosameBunchPileUp = static_cast<bool>(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup));
    }
    isSelColl = isLcFound && isSel8 && isNosameBunchPileUp;
    lcSel(isSelColl);
  }
  PROCESS_SWITCH(HfCorrelatorLcHadronsSelection, processLcSelectionMcRec, "Process Lc Selection McRec", false);

  void processLcSelectionMcGen(aod::McCollision const&,
                               CandidatesLcMcGen const& mcParticles)
  {
    bool isLcFound = true;
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != Pdg::kLambdaCPlus) {
        isLcFound = false;
        continue;
      }
      double yL = RecoDecay::y(particle.pVector(), MassLambdaCPlus);
      if (std::abs(yL) > yCandMax || particle.pt() < ptCandMin) {
        isLcFound = false;
        continue;
      }
      isLcFound = true;
      break;
    }
    lcSel(isLcFound);
  }
  PROCESS_SWITCH(HfCorrelatorLcHadronsSelection, processLcSelectionMcGen, "Process Lc Selection McGen", false);
};

// Lc-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via Mc truth)
struct HfCorrelatorLcHadrons {
  Produces<aod::LcHadronPair> entryLcHadronPair;
  Produces<aod::LcHadronPairTrkPID> entryLcHadronPairTrkPID;
  Produces<aod::LcHadronRecoInfo> entryLcHadronRecoInfo;
  Produces<aod::LcHadronMlInfo> entryLcHadronMlInfo;
  Produces<aod::LcRecoInfo> entryLcCandRecoInfo;
  Produces<aod::LcHadronGenInfo> entryLcHadronGenInfo;
  Produces<aod::LcGenInfo> entryLcCandGenInfo;
  Produces<aod::TrkRecInfoLc> entryTrackRecoInfo;
  Produces<aod::Lc> entryLc;
  Produces<aod::Hadron> entryHadron;
  Produces<aod::LcHadronTrkPID> entryTrkPID;

  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "number of events mixed in ME process"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying Lc efficiency weights"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCAxy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCAz of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "max. cand. pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 50., "max. track pT"};
  Configurable<float> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<float> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<std::vector<double>> binsPtLc{"binsPtLc", std::vector<double>{o2::analysis::hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle"};
  Configurable<std::vector<double>> binsPtEfficiencyLc{"binsPtEfficiencyLc", std::vector<double>{o2::analysis::hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits for efficiency"};
  Configurable<std::vector<double>> efficiencyLc{"efficiencyLc", {1., 1., 1., 1., 1., 1.}, "efficiency values for Lc"};
  Configurable<bool> storeAutoCorrelationFlag{"storeAutoCorrelationFlag", false, "Store flag that indicates if the track is paired to its Lc mother instead of skipping it"};
  Configurable<bool> correlateLcWithLeadingParticle{"correlateLcWithLeadingParticle", false, "Switch for correlation of Lc baryons with leading particle only"};
  Configurable<bool> pidTrkApplied{"pidTrkApplied", false, "Apply PID selection for associated tracks"};
  Configurable<std::vector<int>> trkPIDspecies{"trkPIDspecies", std::vector<int>{o2::track::PID::Proton, o2::track::PID::Pion, o2::track::PID::Kaon}, "Trk sel: Particles species for PID, proton, pion, kaon"};
  Configurable<std::vector<float>> pidTPCMax{"pidTPCMax", std::vector<float>{3., 0., 0.}, "maximum nSigma TPC"};
  Configurable<std::vector<float>> pidTOFMax{"pidTOFMax", std::vector<float>{3., 0., 0.}, "maximum nSigma TOF"};
  Configurable<float> tofPIDThreshold{"tofPIDThreshold", 0.75, "minimum pT after which TOF PID is applicable"};
  Configurable<bool> fillTrkPID{"fillTrkPID", false, "fill PID information for associated tracks"};
  Configurable<bool> forceTOF{"forceTOF", false, "fill PID information for associated tracks"};
  Configurable<bool> calTrkEff{"calTrkEff", false, "fill histograms to calculate efficiency"};
  Configurable<bool> isRecTrkPhyPrimary{"isRecTrkPhyPrimary", true, "Calculate the efficiency of reconstructed primary physical tracks"};
  Configurable<bool> calEffLcEvent{"calEffLcEvent", true, "Calculate the efficiency of Lc candidate"};

  HfHelper hfHelper;
  SliceCache cache;
  Service<o2::framework::O2DatabasePDG> pdg;
  int leadingIndex = 0;
  bool correlationStatus = false;

  // Event Mixing for the Data Mode
  using SelCollisionsWithLc = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::LcSelection>>;
  using SelCollisionsWithLcMc = soa::Filtered<soa::Join<aod::McCollisions, aod::LcSelection, aod::MultsExtraMC>>; // collisionFilter applied
  using CandidatesLcData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  // Event Mixing for the MCRec Mode
  using CandidatesLcMcRec = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi, aod::HfCand3ProngMcRec>>;
  using CandidatesLcMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>; // flagLcFilter applied
  // Event Mixing for the MCGen Mode
  using McCollisionsSel = soa::Filtered<soa::Join<aod::McCollisions, aod::LcSelection>>;
  using McParticlesSel = soa::Filtered<aod::McParticles>;
  // Tracks used in Data and MC
  using TracksData = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;                           // trackFilter applied
  using TracksWithMc = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>; // trackFilter applied
  // Filters for ME
  Filter collisionFilter = aod::hf_selection_lc_collision::lcSel == true;
  Filter lcFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (nabs(aod::track::pt) > ptTrackMin) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);

  // configurable axis definition
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 100000.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis binsZVtx{"binsZVtx", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "Mixing bins - z-vertex"};
  ConfigurableAxis binsMultiplicityMc{"binsMultiplicityMc", {VARIABLE_WIDTH, 0.0f, 20.0f, 50.0f, 500.0f}, "Mixing bins - MC multiplicity"}; // In MCGen multiplicity is defined by counting tracks
  ConfigurableAxis binsBdtScore{"binsBdtScore", {100, 0., 1.}, "Bdt output scores"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 6000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsMassLc{"binsMassLc", {200, 1.98, 2.58}, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};

  BinningType corrBinning{{binsZVtx, binsMultiplicity}, true};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext&)
  {
    AxisSpec axisMassLc = {binsMassLc, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};
    AxisSpec axisEta = {binsEta, "#it{eta}"};
    AxisSpec axisPhi = {binsPhi, "#it{#varphi}"};
    AxisSpec axisPtLc = {(std::vector<double>)binsPtLc, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T} Hadron (GeV/#it{c})"};
    AxisSpec axisMultiplicity = {binsMultiplicity, "Multiplicity"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisPosZ = {binsZVtx, "PosZ"};
    AxisSpec axisBdtScore = {binsBdtScore, "Bdt score"};
    AxisSpec axisPoolBin = {binsPoolBin, "PoolBin"};
    AxisSpec axisRapidity = {100, -2, 2, "Rapidity"};
    AxisSpec axisSign = {2, -1, 1, "Sign"};

    registry.add("hPtCand", "Lc,Hadron candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtProng0", "Lc,Hadron candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtProng1", "Lc,Hadron candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtProng2", "Lc,Hadron candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hSelectionStatusLcToPKPi", "Lc,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}});
    registry.add("hSelectionStatusLcToPiKP", "Lc,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}});
    registry.add("hEta", "Lc,Hadron candidates;candidate #it{#eta};entries", {HistType::kTH1F, {axisEta}});
    registry.add("hPhi", "Lc,Hadron candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {axisPhi}});
    registry.add("hY", "Lc,Hadron candidates;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hCountLcHadronPerEvent", "Lc,Hadron particles - MC gen;Number per event;entries", {HistType::kTH1F, {{21, -0.5, 20.5}}});
    registry.add("hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hMultFT0M", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hZvtx", "z vertex;z vertex;entries", {HistType::kTH1F, {{200, -20., 20.}}});
    registry.add("hLcBin", "Lc selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}});
    registry.add("hTracksBin", "Tracks selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}});
    registry.add("hMassLcVsPt", "Lc candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisMassLc}, {axisPtLc}}});
    registry.add("hMassLcData", "Lc candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{axisMassLc}}});
    registry.add("hLcPoolBin", "Lc candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hTracksPoolBin", "Particles associated pool bin", {HistType::kTH1F, {axisPoolBin}});
    // Histograms for MC Reco analysis
    registry.add("hSelectionStatusLcToPKPiMcRec", "Lc,Hadron candidates - MC reco;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}});
    registry.add("hSelectionStatusLcToPiKPMcRec", "Lc,Hadron candidates - MC reco;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}});
    registry.add("hMcEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}});
    registry.add("hPtProng0McRec", "Lc,Hadron candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtProng1McRec", "Lc,Hadron candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtProng2McRec", "Lc,Hadron candidates - MC reco;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hMassLcMcRec", "Lc candidates;inv. mass (P K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{axisMassLc}}});
    registry.add("hMassLcVsPtMcRec", "Lc candidates - MC Reco", {HistType::kTH2F, {{axisMassLc}, {axisPtLc}}});
    registry.add("hMassLcMcRecSig", "Lc signal candidates - MC reco;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisMassLc}, {axisPtLc}}});
    registry.add("hMassLcMcRecBkg", "Lc background candidates - MC reco;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisMassLc}, {axisPtLc}}});
    registry.add("hPtCandMcRecSig", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtCandMcRecSigPrompt", "Lc,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtCandMcRecSigNonPrompt", "Lc,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtCandMcRecBkg", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtLc}});
    registry.add("hEtaMcRecSig", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcRecSig", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
    registry.add("hYMcRecSig", "Lc,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hEtaMcRecBkg", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcRecBkg", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
    registry.add("hYMcRecBkg", "Lc,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hFakeTracksMcRec", "Fake tracks - MC Rec", {HistType::kTH1F, {axisPtHadron}});
    registry.add("hPtParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtLc}}});
    registry.add("hPtTracksVsSignRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisSign}}});
    registry.add("hPtTracksVsSignRecTrue", "Associated Particle - MC Rec (True)", {HistType::kTH2F, {{axisPtHadron}, {axisSign}}});
    registry.add("hPtTracksVsSignGen", "Associated Particle - MC Gen", {HistType::kTH2F, {{axisPtHadron}, {axisSign}}});
    registry.add("hPtPrimaryParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtLc}}});
    registry.add("hPtVsMultiplicityMcRecPrompt", "Multiplicity FT0M - MC Rec Prompt", {HistType::kTH2F, {{axisPtLc}, {axisMultFT0M}}});
    registry.add("hPtVsMultiplicityMcRecNonPrompt", "Multiplicity FT0M - MC Rec Non Prompt", {HistType::kTH2F, {{axisPtLc}, {axisMultFT0M}}});
    // Histograms for MC Gen analysis
    registry.add("hcountLctriggersMcGen", "Lc trigger particles - MC gen;;N of trigger Lc", {HistType::kTH2F, {{1, -0.5, 0.5}, {axisPtLc}}});
    registry.add("hPtCandMcGen", "Lc,Hadron particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hYMcGen", "Lc,Hadron candidates - MC gen;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hPtCandMcGenPrompt", "Lc,Hadron particles - MC Gen Prompt", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtCandMcGenNonPrompt", "Lc,Hadron particles - MC Gen Non Prompt", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtParticleAssocMcGen", "Associated Particle - MC Gen", {HistType::kTH1F, {axisPtHadron}});
    registry.add("hEtaMcGen", "Lc,Hadron particles - MC Gen", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcGen", "Lc,Hadron particles - MC Gen", {HistType::kTH1F, {axisPhi}});
    registry.add("hMultFT0AMcGen", "Lc,Hadron multiplicity FT0A - MC Gen", {HistType::kTH1F, {axisMultiplicity}});

    corrBinning = {{binsZVtx, binsMultiplicity}, true};
  }

  /// Lc-hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(SelCollisionsWithLc::iterator const& collision,
                   TracksData const& tracks,
                   CandidatesLcData const& candidates, aod::BCsWithTimestamps const&)
  {
    if (candidates.size() == 0) {
      return;
    }

    // find leading particle
    if (correlateLcWithLeadingParticle) {
      leadingIndex = findLeadingParticle(tracks, dcaXYTrackMax.value, dcaZTrackMax.value, etaTrackMax.value);
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int gCollisionId = collision.globalIndex();
    int64_t timeStamp = bc.timestamp();

    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (std::abs(track.eta()) > etaTrackMax || std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
          continue;
        }
        nTracks++;
      }
    }
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < multMin || nTracks > multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    int countLc = 0;
    std::vector<float> outputMl = {-1., -1., -1.};

    for (const auto& candidate : candidates) {
      if (std::abs(hfHelper.yLc(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      double efficiencyWeightLc = 1.;
      if (applyEfficiency) {
        efficiencyWeightLc = 1. / efficiencyLc->at(o2::analysis::findBin(binsPtEfficiencyLc, candidate.pt()));
      }
      auto trackPos1 = candidate.template prong0_as<TracksData>(); // positive daughter (negative for the antiparticles)
      int8_t chargeLc = trackPos1.sign();                          // charge of 1st prong will be the charge of Lc candidate

      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hPtProng2"), candidate.ptProng2());
      registry.fill(HIST("hEta"), candidate.eta());
      registry.fill(HIST("hPhi"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));
      registry.fill(HIST("hY"), hfHelper.yLc(candidate));
      registry.fill(HIST("hLcBin"), poolBin);
      if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
        registry.fill(HIST("hMassLcVsPt"), hfHelper.invMassLcToPKPi(candidate), candidate.pt(), efficiencyWeightLc);
        registry.fill(HIST("hMassLcData"), hfHelper.invMassLcToPKPi(candidate), efficiencyWeightLc);
        registry.fill(HIST("hSelectionStatusLcToPKPi"), candidate.isSelLcToPKPi());
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbLcToPKPi()[classMl->at(iclass)];
        }
        entryLcCandRecoInfo(hfHelper.invMassLcToPKPi(candidate), candidate.pt() * chargeLc, outputMl[0], outputMl[1]); // 0: BkgBDTScore, 1:PromptBDTScore
        entryLc(candidate.phi(), candidate.eta(), candidate.pt(), hfHelper.invMassLcToPKPi(candidate), poolBin, gCollisionId, timeStamp);
      }
      if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
        registry.fill(HIST("hMassLcVsPt"), hfHelper.invMassLcToPiKP(candidate), candidate.pt(), efficiencyWeightLc);
        registry.fill(HIST("hMassLcData"), hfHelper.invMassLcToPiKP(candidate), efficiencyWeightLc);
        registry.fill(HIST("hSelectionStatusLcToPiKP"), candidate.isSelLcToPiKP());
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbLcToPiKP()[classMl->at(iclass)];
        }
        entryLcCandRecoInfo(hfHelper.invMassLcToPiKP(candidate), candidate.pt() * chargeLc, outputMl[0], outputMl[1]); // 0: BkgBDTScore, 1:PromptBDTScore
        entryLc(candidate.phi(), candidate.eta(), candidate.pt(), hfHelper.invMassLcToPiKP(candidate), poolBin, gCollisionId, timeStamp);
      }

      // Lc-Hadron correlation dedicated section
      // if the candidate is a Lc, search for Hadrons and evaluate correlations
      for (const auto& track : tracks) {
        // Remove Lc daughters by checking track indices
        if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
          if (!storeAutoCorrelationFlag) {
            continue;
          }
          correlationStatus = true;
        }
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }
        if (pidTrkApplied) {
          if (!passPIDSelection(track, trkPIDspecies, pidTPCMax, pidTOFMax, tofPIDThreshold, forceTOF))
            continue;
        }
        if (correlateLcWithLeadingParticle) {
          if (track.globalIndex() != leadingIndex) {
            continue;
          }
        }
        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt() * chargeLc,
                            track.pt() * track.sign(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPKPi(candidate), false);
          entryLcHadronGenInfo(false, false, 0);
          entryLcHadronMlInfo(outputMl[0], outputMl[1]);
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
          if (fillTrkPID) {
            entryLcHadronPairTrkPID(track.tpcNSigmaPr(), track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tofNSigmaPr(), track.tofNSigmaKa(), track.tofNSigmaPi());
          }
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt() * chargeLc,
                            track.pt() * track.sign(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPiKP(candidate), false);
          entryLcHadronGenInfo(false, false, 0);
          entryLcHadronMlInfo(outputMl[0], outputMl[1]);
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
          if (fillTrkPID) {
            entryLcHadronPairTrkPID(track.tpcNSigmaPr(), track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tofNSigmaPr(), track.tofNSigmaKa(), track.tofNSigmaPi());
          }
        }
        if (countLc == 0) {
          entryHadron(track.phi(), track.eta(), track.pt() * track.sign(), poolBin, gCollisionId, timeStamp);
          if (fillTrkPID) {
            entryTrkPID(track.tpcNSigmaPr(), track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tofNSigmaPr(), track.tofNSigmaKa(), track.tofNSigmaPi());
          }
          registry.fill(HIST("hTracksBin"), poolBin);
        }
      } // Hadron Tracks loop
      countLc++;
    } // end outer Lc loop
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), collision.multFT0M());
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processData, "Process data", true);

  /// Lc-Hadron correlation process starts for McRec
  void processMcRec(SelCollisionsWithLc::iterator const& collision,
                    TracksWithMc const& tracks,
                    CandidatesLcMcRec const& candidates,
                    aod::McParticles const& mcParticles)
  {
    if (candidates.size() == 0) {
      return;
    }

    // find leading particle
    if (correlateLcWithLeadingParticle) {
      leadingIndex = findLeadingParticle(tracks, dcaXYTrackMax.value, dcaZTrackMax.value, etaTrackMax.value);
    }

    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (std::abs(track.eta()) >= etaTrackMax || std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
          continue;
        }
        nTracks++;
      }
    }
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < multMin || nTracks > multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    float multiplicityFT0M = collision.multFT0M();
    // Mc reco level
    bool isLcPrompt = false;
    bool isLcNonPrompt = false;
    bool isLcSignal = false;
    int countLc = 1;
    for (const auto& candidate : candidates) {
      // check decay channel flag for candidate
      if (std::abs(hfHelper.yLc(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      double efficiencyWeightLc = 1.;
      if (applyEfficiency) {
        efficiencyWeightLc = 1. / efficiencyLc->at(o2::analysis::findBin(binsPtEfficiencyLc, candidate.pt()));
      }
      auto trackPos1 = candidate.template prong0_as<TracksWithMc>(); // positive daughter (negative for the antiparticles)
      int8_t chargeLc = trackPos1.sign();                            // charge of 1st prong will be the charge of Lc candidate
      isLcSignal = TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_3prong::DecayType::LcToPKPi);
      isLcPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
      isLcNonPrompt = candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
      std::vector<float> outputMl = {-1., -1., -1.};
      if (isLcSignal) {
        registry.fill(HIST("hPtProng0McRec"), candidate.ptProng0());
        registry.fill(HIST("hPtProng1McRec"), candidate.ptProng1());
        registry.fill(HIST("hPtProng2McRec"), candidate.ptProng2());
        registry.fill(HIST("hPtCandMcRecSig"), candidate.pt());
        registry.fill(HIST("hEtaMcRecSig"), candidate.eta());
        registry.fill(HIST("hPhiMcRecSig"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));
        registry.fill(HIST("hYMcRecSig"), hfHelper.yLc(candidate));
        // LcToPKPi and LcToPiKP division
        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbLcToPKPi()[classMl->at(iclass)];
          }
          // prompt and non-prompt division
          if (isLcPrompt) {
            registry.fill(HIST("hPtCandMcRecSigPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), candidate.pt(), multiplicityFT0M);
          } else if (isLcNonPrompt) {
            registry.fill(HIST("hPtCandMcRecSigNonPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), candidate.pt(), multiplicityFT0M);
          }
          registry.fill(HIST("hMassLcMcRec"), hfHelper.invMassLcToPKPi(candidate), efficiencyWeightLc);
          registry.fill(HIST("hMassLcMcRecSig"), hfHelper.invMassLcToPKPi(candidate), candidate.pt(), efficiencyWeightLc);
          registry.fill(HIST("hMassLcVsPtMcRec"), hfHelper.invMassLcToPKPi(candidate), candidate.pt(), efficiencyWeightLc);
          registry.fill(HIST("hSelectionStatusLcToPKPiMcRec"), candidate.isSelLcToPKPi());
          entryLcCandRecoInfo(hfHelper.invMassLcToPKPi(candidate), candidate.pt() * chargeLc, outputMl[0], outputMl[1]); // 0: BkgBDTScore, 1:PromptBDTScore
          entryLcCandGenInfo(isLcPrompt);
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbLcToPiKP()[classMl->at(iclass)];
          }
          if (isLcPrompt) {
            registry.fill(HIST("hPtCandMcRecSigPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), candidate.pt(), multiplicityFT0M);
          } else if (isLcNonPrompt) {
            registry.fill(HIST("hPtCandMcRecSigNonPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), candidate.pt(), multiplicityFT0M);
          }
          registry.fill(HIST("hMassLcMcRec"), hfHelper.invMassLcToPiKP(candidate), efficiencyWeightLc);
          registry.fill(HIST("hMassLcMcRecSig"), hfHelper.invMassLcToPiKP(candidate), candidate.pt(), efficiencyWeightLc);
          registry.fill(HIST("hMassLcVsPtMcRec"), hfHelper.invMassLcToPiKP(candidate), candidate.pt(), efficiencyWeightLc);
          registry.fill(HIST("hSelectionStatusLcToPiKPMcRec"), candidate.isSelLcToPiKP());
          entryLcCandRecoInfo(hfHelper.invMassLcToPiKP(candidate), candidate.pt() * chargeLc, outputMl[0], outputMl[1]); // 0: BkgBDTScore, 1:PromptBDTScore
          entryLcCandGenInfo(isLcPrompt);
        }
      } else {
        registry.fill(HIST("hPtCandMcRecBkg"), candidate.pt());
        registry.fill(HIST("hEtaMcRecBkg"), candidate.eta());
        registry.fill(HIST("hPhiMcRecBkg"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));
        registry.fill(HIST("hYMcRecBkg"), hfHelper.yLc(candidate));
        // LcToPKPi and LcToPiKP division
        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbLcToPKPi()[classMl->at(iclass)];
          }
          registry.fill(HIST("hMassLcMcRec"), hfHelper.invMassLcToPKPi(candidate), efficiencyWeightLc);
          registry.fill(HIST("hMassLcMcRecBkg"), hfHelper.invMassLcToPKPi(candidate), candidate.pt(), efficiencyWeightLc);
          registry.fill(HIST("hMassLcVsPtMcRec"), hfHelper.invMassLcToPKPi(candidate), candidate.pt(), efficiencyWeightLc);
          registry.fill(HIST("hSelectionStatusLcToPKPiMcRec"), candidate.isSelLcToPKPi());
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbLcToPiKP()[classMl->at(iclass)];
          }
          registry.fill(HIST("hMassLcMcRec"), hfHelper.invMassLcToPiKP(candidate), efficiencyWeightLc);
          registry.fill(HIST("hMassLcMcRecBkg"), hfHelper.invMassLcToPiKP(candidate), candidate.pt(), efficiencyWeightLc);
          registry.fill(HIST("hMassLcVsPtMcRec"), hfHelper.invMassLcToPiKP(candidate), candidate.pt(), efficiencyWeightLc);
          registry.fill(HIST("hSelectionStatusLcToPiKPMcRec"), candidate.isSelLcToPiKP());
        }
      }
      registry.fill(HIST("hLcBin"), poolBin);

      if (calTrkEff && !isLcSignal && calEffLcEvent)
        continue;

      if (calTrkEff && countLc == 1) {
        // genrated tracks
        for (const auto& track : mcParticles) {
          if (std::abs(track.eta()) > etaTrackMax || track.pt() < ptTrackMin || track.pt() > ptTrackMax) {
            continue;
          }
          if ((std::abs(track.pdgCode()) != kElectron) && (std::abs(track.pdgCode()) != kMuonMinus) && (std::abs(track.pdgCode()) != kPiPlus) && (std::abs(track.pdgCode()) != kKPlus) && (std::abs(track.pdgCode()) != kProton)) {
            continue;
          }

          if (pidTrkApplied && (std::abs(track.pdgCode()) != kProton))
            continue; // proton PID

          if (!track.isPhysicalPrimary()) {
            continue;
          }

          int8_t chargeTrack = pdg->GetParticle(track.pdgCode())->Charge(); // Retrieve charge
          registry.fill(HIST("hPtTracksVsSignGen"), track.pt(), chargeTrack);
        }
      }

      // Lc-Hadron correlation dedicated section
      // if the candidate is selected as Lc, search for Hadron ad evaluate correlations
      for (const auto& track : tracks) {
        bool isPhysicalPrimary = false;
        int trackOrigin = -1;
        // apply track selection
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }
        if (pidTrkApplied) {
          if (!passPIDSelection(track, trkPIDspecies, pidTPCMax, pidTOFMax, tofPIDThreshold, forceTOF))
            continue;
        }

        if (calTrkEff && countLc == 1 && track.has_mcParticle()) {
          auto mcParticle = track.template mcParticle_as<aod::McParticles>();
          if (!mcParticle.isPhysicalPrimary() && isRecTrkPhyPrimary)
            continue;

          registry.fill(HIST("hPtTracksVsSignRec"), track.pt(), track.sign());
          if (std::abs(mcParticle.pdgCode()) == kProton)
            registry.fill(HIST("hPtTracksVsSignRecTrue"), track.pt(), track.sign());
        }

        // Removing Lc daughters by checking track indices
        if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
          if (!storeAutoCorrelationFlag) {
            continue;
          }
          correlationStatus = true;
        }

        if (correlateLcWithLeadingParticle) {
          if (track.globalIndex() != leadingIndex) {
            continue;
          }
        }

        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt() * chargeLc,
                            track.pt() * track.sign(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPKPi(candidate), isLcSignal);
          if (fillTrkPID) {
            entryLcHadronPairTrkPID(track.tpcNSigmaPr(), track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tofNSigmaPr(), track.tofNSigmaKa(), track.tofNSigmaPi());
          }
          entryLcHadronMlInfo(outputMl[0], outputMl[1]);
          if (track.has_mcParticle()) {
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            isPhysicalPrimary = mcParticle.isPhysicalPrimary();
            trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
            entryLcHadronGenInfo(isLcPrompt, isPhysicalPrimary, trackOrigin);
          } else {
            entryLcHadronGenInfo(isLcPrompt, false, 0);
            registry.fill(HIST("hFakeTracksMcRec"), track.pt());
          }

          // for secondary particle fraction estimation
          registry.fill(HIST("hPtParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
          if (isPhysicalPrimary) {
            registry.fill(HIST("hPtPrimaryParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
          }
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt() * chargeLc,
                            track.pt() * track.sign(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPiKP(candidate), isLcSignal);
          if (fillTrkPID) {
            entryLcHadronPairTrkPID(track.tpcNSigmaPr(), track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tofNSigmaPr(), track.tofNSigmaKa(), track.tofNSigmaPi());
          }
          entryLcHadronMlInfo(outputMl[0], outputMl[1]);
          if (track.has_mcParticle()) {
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            isPhysicalPrimary = mcParticle.isPhysicalPrimary();
            trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
            entryLcHadronGenInfo(isLcPrompt, isPhysicalPrimary, trackOrigin);
          } else {
            entryLcHadronGenInfo(isLcPrompt, false, 0);
            registry.fill(HIST("hFakeTracksMcRec"), track.pt());
          }
          // for secondary particle fraction estimation
          registry.fill(HIST("hPtParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
          if (isPhysicalPrimary) {
            registry.fill(HIST("hPtPrimaryParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
          }
          entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
        }
      } // end inner loop (Tracks)
      countLc++;
    } // end outer Lc loop
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), collision.multFT0M());
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processMcRec, "Process Mc Reco mode", false);

  /// Lc-Hadron correlation pair builder - for Mc Gen-level analysis
  void processMcGen(SelCollisionsWithLcMc::iterator const& mcCollision,
                    CandidatesLcMcGen const& mcParticles)
  {
    int counterLcHadron = 0;
    registry.fill(HIST("hMcEvtCount"), 0);

    BinningTypeMcGen corrBinningMcGen{{binsZVtx, binsMultiplicityMc}, true};
    int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), mcCollision.multMCFT0A()));
    registry.fill(HIST("hMultFT0AMcGen"), mcCollision.multMCFT0A());

    bool isLcPrompt = false;
    bool isLcNonPrompt = false;

    // find leading particle
    if (correlateLcWithLeadingParticle) {
      leadingIndex = findLeadingParticleMcGen(mcParticles, etaTrackMax.value, ptTrackMin.value);
    }

    // Mc Gen level
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != Pdg::kLambdaCPlus) {
        continue;
      }
      if (!TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        continue;
      }
      double yL = RecoDecay::y(particle.pVector(), MassLambdaCPlus);
      if (std::abs(yL) > yCandMax || particle.pt() < ptCandMin) {
        continue;
      }
      registry.fill(HIST("hLcBin"), poolBin);
      registry.fill(HIST("hPtCandMcGen"), particle.pt());
      registry.fill(HIST("hEtaMcGen"), particle.eta());
      registry.fill(HIST("hPhiMcGen"), RecoDecay::constrainAngle(particle.phi(), -PIHalf));
      registry.fill(HIST("hYMcGen"), yL);

      isLcPrompt = particle.originMcGen() == RecoDecay::OriginType::Prompt;
      isLcNonPrompt = particle.originMcGen() == RecoDecay::OriginType::NonPrompt;
      if (isLcPrompt) {
        registry.fill(HIST("hPtCandMcGenPrompt"), particle.pt());
      } else if (isLcNonPrompt) {
        registry.fill(HIST("hPtCandMcGenNonPrompt"), particle.pt());
      }

      // prompt and non-prompt division
      std::vector<int> listDaughters{};
      std::array<int, 3> arrDaughLcPDG = {kProton, -kKPlus, kPiPlus};
      std::array<int, 3> prongsId;
      listDaughters.clear();
      RecoDecay::getDaughters(particle, &listDaughters, arrDaughLcPDG, 2);
      int counterDaughters = 0;
      if (listDaughters.size() == 3) {
        for (const auto& dauIdx : listDaughters) {
          auto daughI = mcParticles.rawIteratorAt(dauIdx - mcParticles.offset());
          counterDaughters += 1;
          prongsId[counterDaughters - 1] = daughI.globalIndex();
        }
      }
      counterLcHadron++;
      // Lc Hadron correlation dedicated section
      // if it's a Lc particle, search for Hadron and evalutate correlations
      registry.fill(HIST("hcountLctriggersMcGen"), 0, particle.pt()); // to count trigger Lc for normalisation
      for (const auto& particleAssoc : mcParticles) {
        if (std::abs(particleAssoc.eta()) > etaTrackMax || particleAssoc.pt() < ptTrackMin || particleAssoc.pt() > ptTrackMax) {
          continue;
        }
        if (particleAssoc.globalIndex() == prongsId[0] || particleAssoc.globalIndex() == prongsId[1] || particleAssoc.globalIndex() == prongsId[2]) {
          if (!storeAutoCorrelationFlag) {
            continue;
          }
          correlationStatus = true;
        }

        if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particle.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
          continue;
        }

        if (pidTrkApplied && (std::abs(particleAssoc.pdgCode()) != kProton))
          continue; // proton PID

        if (!particleAssoc.isPhysicalPrimary()) {
          continue;
        }

        if (correlateLcWithLeadingParticle) {
          if (particleAssoc.globalIndex() != leadingIndex) {
            continue;
          }
        }

        int8_t chargeLc = pdg->GetParticle(particle.pdgCode())->Charge();         // Retrieve charge
        int8_t chargeAssoc = pdg->GetParticle(particleAssoc.pdgCode())->Charge(); // Retrieve charge

        int trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, particleAssoc, true);
        registry.fill(HIST("hPtParticleAssocMcGen"), particleAssoc.pt());
        entryLcHadronPair(getDeltaPhi(particleAssoc.phi(), particle.phi()),
                          particleAssoc.eta() - particle.eta(),
                          particle.pt() * chargeLc,
                          particleAssoc.pt() * chargeAssoc,
                          poolBin,
                          correlationStatus);
        entryLcHadronRecoInfo(MassLambdaCPlus, true);
        entryLcHadronGenInfo(isLcPrompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
      } // end inner loop
    } // end outer loop
    registry.fill(HIST("hCountLcHadronPerEvent"), counterLcHadron);
    registry.fill(HIST("hZvtx"), mcCollision.posZ());
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processMcGen, "Process Mc Gen mode", false);

  void processDataMixedEvent(SelCollisionsWithLc const& collisions,
                             CandidatesLcData const& candidates,
                             TracksData const& tracks)
  {
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelCollisionsWithLc, CandidatesLcData, TracksData, BinningType> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      for (const auto& [trigLc, assocParticle] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!assocParticle.isGlobalTrackWoDCA() || std::abs(hfHelper.yLc(trigLc)) > yCandMax) {
          continue;
        }

        if (pidTrkApplied) {
          if (!passPIDSelection(assocParticle, trkPIDspecies, pidTPCMax, pidTOFMax, tofPIDThreshold, forceTOF))
            continue;
        }

        auto trackPos1 = trigLc.template prong0_as<TracksData>(); // positive daughter (negative for the antiparticles)
        int8_t chargeLc = trackPos1.sign();                       // charge of 1st prong will be the charge of Lc candidate

        std::vector<float> outputMl = {-1., -1., -1.};
        // LcToPKPi and LcToPiKP division
        if (trigLc.isSelLcToPKPi() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(assocParticle.phi(), trigLc.phi()),
                            assocParticle.eta() - trigLc.eta(),
                            trigLc.pt() * chargeLc,
                            assocParticle.pt() * assocParticle.sign(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPKPi(trigLc), false);
          entryLcHadronGenInfo(false, false, 0);
          if (fillTrkPID) {
            entryLcHadronPairTrkPID(assocParticle.tpcNSigmaPr(), assocParticle.tpcNSigmaKa(), assocParticle.tpcNSigmaPi(), assocParticle.tofNSigmaPr(), assocParticle.tofNSigmaKa(), assocParticle.tofNSigmaPi());
          }
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = trigLc.mlProbLcToPKPi()[classMl->at(iclass)];
          }
          entryLcHadronMlInfo(outputMl[0], outputMl[1]);
          entryTrackRecoInfo(assocParticle.dcaXY(), assocParticle.dcaZ(), assocParticle.tpcNClsCrossedRows());
        }
        if (trigLc.isSelLcToPiKP() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(assocParticle.phi(), trigLc.phi()),
                            assocParticle.eta() - trigLc.eta(),
                            trigLc.pt() * chargeLc,
                            assocParticle.pt() * assocParticle.sign(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPiKP(trigLc), false);
          entryLcHadronGenInfo(false, false, 0);
          if (fillTrkPID) {
            entryLcHadronPairTrkPID(assocParticle.tpcNSigmaPr(), assocParticle.tpcNSigmaKa(), assocParticle.tpcNSigmaPi(), assocParticle.tofNSigmaPr(), assocParticle.tofNSigmaKa(), assocParticle.tofNSigmaPi());
          }
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = trigLc.mlProbLcToPiKP()[classMl->at(iclass)];
          }
          entryLcHadronMlInfo(outputMl[0], outputMl[1]);
          entryTrackRecoInfo(assocParticle.dcaXY(), assocParticle.dcaZ(), assocParticle.tpcNClsCrossedRows());
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processDataMixedEvent, "Process Mixed Event Data", false);

  void processMcRecMixedEvent(SelCollisionsWithLc const& collisions,
                              CandidatesLcMcRec const& candidates,
                              TracksWithMc const& tracks,
                              aod::McParticles const& mcParticles)
  {
    BinningType corrBinning{{binsZVtx, binsMultiplicityMc}, true};
    for (const auto& candidate : candidates) {
      if (std::abs(hfHelper.yLc(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      // Lc flag
      bool isLcSignal = TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_3prong::DecayType::LcToPKPi);
      // prompt and non-prompt division
      bool isLcPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
      bool isLcNonPrompt = candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
      if (isLcSignal) {
        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          if (isLcPrompt) {
            registry.fill(HIST("hPtCandMcRecSigPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), candidate.pt(), 0);
          } else if (isLcNonPrompt) {
            registry.fill(HIST("hPtCandMcRecSigNonPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), candidate.pt(), 0);
          }
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          if (isLcPrompt) {
            registry.fill(HIST("hPtCandMcRecSigPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), candidate.pt(), 0);
          } else if (isLcNonPrompt) {
            registry.fill(HIST("hPtCandMcRecSigNonPrompt"), candidate.pt());
            registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), candidate.pt(), 0);
          }
        }
      } else {
        registry.fill(HIST("hPtCandMcRecBkg"), candidate.pt());
        registry.fill(HIST("hEtaMcRecBkg"), candidate.eta());
        registry.fill(HIST("hPhiMcRecBkg"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));
      }
    }
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelCollisionsWithLc, CandidatesLcMcRec, TracksWithMc, BinningType> pairMcRec{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairMcRec) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      int poolBinLc = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M()));
      registry.fill(HIST("hMultFT0M"), c1.multFT0M());
      registry.fill(HIST("hZvtx"), c1.posZ());
      registry.fill(HIST("hTracksPoolBin"), poolBin);
      registry.fill(HIST("hLcPoolBin"), poolBinLc);
      for (const auto& [candidate, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (std::abs(hfHelper.yLc(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
          continue;
        }
        if (!pAssoc.isGlobalTrackWoDCA()) {
          continue;
        }
        std::vector<float> outputMl = {-1., -1., -1.};
        bool isPhysicalPrimary = false;
        int trackOrigin = -1;
        bool isLcSignal = std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi;
        bool isLcPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
        if (pidTrkApplied) {
          if (!passPIDSelection(pAssoc, trkPIDspecies, pidTPCMax, pidTOFMax, tofPIDThreshold, forceTOF))
            continue;
        }
        auto trackPos1 = candidate.template prong0_as<TracksWithMc>(); // positive daughter (negative for the antiparticles)
        int8_t chargeLc = trackPos1.sign();                            // charge of 1st prong will be the charge of Lc candidate

        if (pAssoc.has_mcParticle()) {
          auto mcParticle = pAssoc.template mcParticle_as<aod::McParticles>();
          isPhysicalPrimary = mcParticle.isPhysicalPrimary();
          trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
        } else {
          registry.fill(HIST("hFakeTracksMcRec"), pAssoc.pt());
        }
        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(pAssoc.phi(), candidate.phi()),
                            pAssoc.eta() - candidate.eta(),
                            candidate.pt() * chargeLc,
                            pAssoc.pt() * pAssoc.sign(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPKPi(candidate), isLcSignal);
          entryLcHadronGenInfo(isLcPrompt, isPhysicalPrimary, trackOrigin);
          if (fillTrkPID) {
            entryLcHadronPairTrkPID(pAssoc.tpcNSigmaPr(), pAssoc.tpcNSigmaKa(), pAssoc.tpcNSigmaPi(), pAssoc.tofNSigmaPr(), pAssoc.tofNSigmaKa(), pAssoc.tofNSigmaPi());
          }
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbLcToPKPi()[classMl->at(iclass)];
          }
          entryLcHadronMlInfo(outputMl[0], outputMl[1]);
          entryTrackRecoInfo(pAssoc.dcaXY(), pAssoc.dcaZ(), pAssoc.tpcNClsCrossedRows());
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(pAssoc.phi(), candidate.phi()),
                            pAssoc.eta() - candidate.eta(),
                            candidate.pt() * chargeLc,
                            pAssoc.pt() * pAssoc.sign(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPiKP(candidate), isLcSignal);
          entryLcHadronGenInfo(isLcPrompt, isPhysicalPrimary, trackOrigin);
          if (fillTrkPID) {
            entryLcHadronPairTrkPID(pAssoc.tpcNSigmaPr(), pAssoc.tpcNSigmaKa(), pAssoc.tpcNSigmaPi(), pAssoc.tofNSigmaPr(), pAssoc.tofNSigmaKa(), pAssoc.tofNSigmaPi());
          }
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbLcToPiKP()[classMl->at(iclass)];
          }
          entryLcHadronMlInfo(outputMl[0], outputMl[1]);
          entryTrackRecoInfo(pAssoc.dcaXY(), pAssoc.dcaZ(), pAssoc.tpcNClsCrossedRows());
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processMcRecMixedEvent, "Process Mixed Event McRec", false);

  void processMcGenMixedEvent(SelCollisionsWithLcMc const& collisions,
                              CandidatesLcMcGen const& mcParticles)
  {
    BinningTypeMcGen corrBinningMcGen{{binsZVtx, binsMultiplicityMc}, true};
    auto tracksTuple = std::make_tuple(mcParticles, mcParticles);
    Pair<SelCollisionsWithLcMc, CandidatesLcMcGen, CandidatesLcMcGen, BinningTypeMcGen> pairMcGen{corrBinningMcGen, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcGen) {
      int poolBin = corrBinningMcGen.getBin(std::make_tuple(c1.posZ(), c1.multMCFT0A()));
      for (const auto& [candidate, particleAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (std::abs(candidate.pdgCode()) != Pdg::kLambdaCPlus) {
          continue;
        }
        double yL = RecoDecay::y(candidate.pVector(), MassLambdaCPlus);
        if (std::abs(yL) > yCandGenMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
          continue;
        }
        if (std::abs(particleAssoc.eta()) > etaTrackMax || particleAssoc.pt() < ptTrackMin || particleAssoc.pt() > ptTrackMax) {
          continue;
        }
        if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particleAssoc.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
          continue;
        }
        if (!particleAssoc.isPhysicalPrimary()) {
          continue;
        }
        if (pidTrkApplied && (std::abs(particleAssoc.pdgCode()) != kProton)) {
          continue; // proton PID
        }
        int8_t chargeLc = pdg->GetParticle(candidate.pdgCode())->Charge();        // Retrieve charge
        int8_t chargeAssoc = pdg->GetParticle(particleAssoc.pdgCode())->Charge(); // Retrieve charge

        int trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, particleAssoc, true);
        bool isLcPrompt = candidate.originMcGen() == RecoDecay::OriginType::Prompt;
        entryLcHadronPair(getDeltaPhi(particleAssoc.phi(), candidate.phi()),
                          particleAssoc.eta() - candidate.eta(),
                          candidate.pt() * chargeLc,
                          particleAssoc.pt() * chargeAssoc,
                          poolBin,
                          correlationStatus);
        entryLcHadronRecoInfo(MassLambdaCPlus, true);
        entryLcHadronGenInfo(isLcPrompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processMcGenMixedEvent, "Process Mixed Event McGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorLcHadronsSelection>(cfgc),
                      adaptAnalysisTask<HfCorrelatorLcHadrons>(cfgc)};
}
