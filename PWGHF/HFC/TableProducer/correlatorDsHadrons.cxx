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

/// \file correlatorDsHadrons.cxx
/// \brief Ds-Hadrons correlator task - data-like, MC-reco and MC-Gen analyses
/// \author Grazia Luparello <grazia.luparello@cern.ch>
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
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

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
double getDeltaPhi(double phiHadron, double phiD)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -PIHalf);
}

// binning type
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;

/// Code to select collisions with at least one Ds meson
struct HfCorrelatorDsHadronsSelCollision {
  Produces<aod::DmesonSelection> collisionsWithSelDs;

  Configurable<bool> useSel8{"useSel8", true, "Flag for applying sel8 for collision selection"};
  Configurable<bool> doSelDsCollision{"doSelDsCollision", true, "Select collisions with at least one Ds"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  SliceCache cache;
  HfHelper hfHelper;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Filter dsFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs);

  /// Code to select collisions with at least one Ds meson - for real data and data-like analysis
  void processDsSelCollisionsData(SelCollisions::iterator const& collision,
                                  CandDsData const& candidates)
  {
    bool isDsFound = false;
    bool isSel8 = false;
    if (doSelDsCollision) {
      for (const auto& candidate : candidates) {
        if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          continue;
        }
        isDsFound = true;
        break;
      }
    } else {
      isDsFound = true;
    }
    if (useSel8) {
      isSel8 = collision.sel8();
      isDsFound = isDsFound && isSel8;
    }
    collisionsWithSelDs(isDsFound);
  }
  PROCESS_SWITCH(HfCorrelatorDsHadronsSelCollision, processDsSelCollisionsData, "Process Ds Collision Selection Data", true);

  /// Code to select collisions with at least one Ds meson - for MC reco-level analysis
  void processDsSelCollisionsMcRec(SelCollisions::iterator const& collision,
                                   CandDsMcReco const& candidates)
  {
    bool isDsFound = false;
    bool isSel8 = false;
    if (doSelDsCollision) {
      for (const auto& candidate : candidates) {
        if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          continue;
        }
        isDsFound = true;
        break;
      }
    } else {
      isDsFound = true;
    }
    if (useSel8) {
      isSel8 = collision.sel8();
      isDsFound = isDsFound && isSel8;
    }
    collisionsWithSelDs(isDsFound);
  }
  PROCESS_SWITCH(HfCorrelatorDsHadronsSelCollision, processDsSelCollisionsMcRec, "Process Ds Collision Selection MCRec", false);

  /// Code to select collisions with at least one Ds meson - for MC gen-level analysis
  void processDsSelCollisionsMcGen(aod::McCollision const&,
                                   CandDsMcGen const& mcParticles)
  {
    bool isDsFound = false;
    if (doSelDsCollision) {
      for (const auto& particle : mcParticles) {
        if (std::abs(particle.pdgCode()) != Pdg::kDS) {
          continue;
        }
        double yD = RecoDecay::y(particle.pVector(), MassDS);
        if (std::abs(yD) > yCandMax || particle.pt() < ptCandMin) {
          continue;
        }
        isDsFound = true;
        break;
      }
    } else {
      isDsFound = true;
    }
    collisionsWithSelDs(isDsFound);
  }
  PROCESS_SWITCH(HfCorrelatorDsHadronsSelCollision, processDsSelCollisionsMcGen, "Process Ds Collision Selection MCGen", false);
};

/// Ds-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorDsHadrons {
  Produces<aod::DsHadronPair> entryDsHadronPair;
  Produces<aod::DsHadronRecoInfo> entryDsHadronRecoInfo;
  Produces<aod::DsHadronGenInfo> entryDsHadronGenInfo;
  Produces<aod::DsHadronMlInfo> entryDsHadronMlInfo;
  Produces<aod::DsCandRecoInfo> entryDsCandRecoInfo;

  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<bool> useSel8ForTrackEff{"useSel8ForTrackEff", true, "Flag for applying sel8 for collision selection"};
  Configurable<bool> applyEfficiency{"applyEfficiency", true, "Flag for applying D-meson efficiency weights"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCA_xy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCA_z of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "max. cand pT"};
  Configurable<float> ptDaughterMin{"ptDaughterMin", 0.1, "min. daughter pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 50., "max. track pT"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle"};
  Configurable<std::vector<double>> binsPtEfficiencyD{"binsPtEfficiencyD", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", {1., 1., 1., 1., 1., 1.}, "efficiency values for Ds meson"};
  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0}, "z vertex position pools"};
  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 900., 1800., 6000.}, "event multiplicity pools (FT0M)"};
  ConfigurableAxis binsMassD{"binsMassD", {200, 1.7, 2.25}, "inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {200, 0., 800.}, "Multiplicity"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 6000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "primary vertex z coordinate"};
  ConfigurableAxis binsBdtScore{"binsBdtScore", {100, 0., 1.}, "Bdt output scores"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};

  HfHelper hfHelper;
  SliceCache cache;

  enum CandidateStep { kCandidateStepMcGenAll = 0,
                       kCandidateStepMcGenDsToKKPi,
                       kCandidateStepMcCandInAcceptance,
                       kCandidateStepMcDaughtersInAcceptance,
                       kCandidateStepMcReco,
                       kCandidateStepMcRecoInAcceptance,
                       kCandidateNSteps };

  enum AssocTrackStep { kAssocTrackStepMcGen = 0,
                        kAssocTrackStepMcGenInAcceptance,
                        kAssocTrackStepRecoAll,
                        kAssocTrackStepRecoMcMatch,
                        kAssocTrackStepRecoPrimaries,
                        kAssocTrackStepFake,
                        kAssocTrackNSteps };

  using SelCollisionsWithDs = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::DmesonSelection>>;            // collisionFilter applied
  using SelCollisionsWithDsMc = soa::Filtered<soa::Join<aod::McCollisions, aod::DmesonSelection>>;                                 // collisionFilter applied
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;                           // flagDsFilter applied
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi, aod::HfCand3ProngMcRec>>; // flagDsFilter applied
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;                                                         // flagDsFilter applied
  using MyTracksData = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection>>;                                             // trackFilter applied
  using TracksWithMc = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, o2::aod::McTrackLabels>>;                     // trackFilter applied

  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter flagDsFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs);
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    AxisSpec axisMassD = {binsMassD, "inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
    AxisSpec axisEta = {binsEta, "#it{#eta}"};
    AxisSpec axisPhi = {binsPhi, "#it{#varphi}"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtD, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T} Hadron (GeV/#it{c})"};
    AxisSpec axisMultiplicity = {binsMultiplicity, "Multiplicity"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    AxisSpec axisBdtScore = {binsBdtScore, "Bdt score"};
    AxisSpec axisPoolBin = {binsPoolBin, "PoolBin"};
    AxisSpec axisStatus = {15, 0.5, 15.5, "Selection status"};

    // Histograms for data analysis
    registry.add("hPtCand", "Ds candidates pt", {HistType::kTH1F, {axisPtD}});
    registry.add("hSelectionStatusDsToKKPi", "Ds candidates selection", {HistType::kTH1F, {axisStatus}});
    registry.add("hSelectionStatusDsToPiKK", "Ds candidates selection", {HistType::kTH1F, {axisStatus}});
    registry.add("hCountSelectionStatusDsToKKPiAndToPiKK", "Ds candidates selection", {HistType::kTH1F, {{1, -0.5, 0.5}}});
    registry.add("hEta", "Ds candidates eta", {HistType::kTH1F, {axisEta}});
    registry.add("hEtaVsPtCand", "Ds candidates etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtD}}});
    registry.add("hEtaVsPtPartAssoc", "Particles associated etaVsPt", {HistType::kTH2F, {{axisEta}, {axisPtHadron}}});
    registry.add("hPhi", "Ds candidates phi", {HistType::kTH1F, {axisPhi}});
    registry.add("hPhiVsPtCand", "Ds candidates phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtD}}});
    registry.add("hPhiVsPtPartAssoc", "Particles associated phiVsPt", {HistType::kTH2F, {{axisPhi}, {axisPtHadron}}});
    registry.add("hMultiplicity", "Multiplicity", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("hMultFT0M", "Multiplicity FT0M", {HistType::kTH1F, {axisMultFT0M}});
    registry.add("hZVtx", "z vertex", {HistType::kTH1F, {axisPosZ}});
    registry.add("hMassDsVsPt", "Ds candidates massVsPt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassDsData", "Ds candidates mass", {HistType::kTH1F, {axisMassD}});
    registry.add("hCollisionPoolBin", "Ds candidates collision pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hDsPoolBin", "Ds candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hTracksPoolBin", "Particles associated pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hCorrelSystematics", "Ds-h correlations systematic error evaluation", {HistType::kTHnSparseD, {{axisPhi}, {axisEta}, {axisPtD}, {axisPtHadron}, {axisMassD}, {axisBdtScore}, {axisBdtScore}}});
    // Histograms for MC Reco analysis
    registry.add("hPtCandMcRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcRecSigPrompt", "Ds,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcRecSigNonPrompt", "Ds,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtParticleAssocMcRec", "Associated Particle - MC Rec", {HistType::kTH1F, {axisPtHadron}});
    registry.add("hPtParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtD}}});
    registry.add("hPtPrimaryParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtD}}});
    registry.add("hEtaMcRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
    registry.add("hEtaMcRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
    registry.add("hPtVsMultiplicityMcRecPrompt", "Multiplicity FT0M - MC Rec Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultFT0M}}});
    registry.add("hPtVsMultiplicityMcRecNonPrompt", "Multiplicity FT0M - MC Rec Non Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultFT0M}}});
    registry.add("hMassDsMcRec", "Ds candidates", {HistType::kTH1F, {axisMassD}});
    registry.add("hMassDsVsPtMcRec", "Ds signal candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassDsMcRecSig", "Ds signal candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassDsMcRecBkg", "Ds background candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hFakeTracksMcRec", "Fake tracks - MC Rec", {HistType::kTH1F, {axisPtHadron}});
    // Histograms for MC Gen analysis
    registry.add("hPtCandMcGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcGenPrompt", "Ds,Hadron particles - MC Gen Prompt", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcGenNonPrompt", "Ds,Hadron particles - MC Gen Non Prompt", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtParticleAssocMcGen", "Associated Particle - MC Gen", {HistType::kTH1F, {axisPtHadron}});
    registry.add("hEtaMcGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisPhi}});
    // Histograms for efficiencies
    auto hCandidates = registry.add<StepTHn>("hCandidates", "Candidate count at different steps", {HistType::kStepTHnF, {axisPtD, axisMultFT0M, {RecoDecay::OriginType::NonPrompt + 1, +RecoDecay::OriginType::None - 0.5, +RecoDecay::OriginType::NonPrompt + 0.5}}, kCandidateNSteps});
    hCandidates->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hCandidates->GetAxis(1)->SetTitle("multiplicity");
    hCandidates->GetAxis(2)->SetTitle("Charm hadron origin");
    auto hAssocTracks = registry.add<StepTHn>("hAssocTracks", "Associated tracks at different steps", {HistType::kStepTHnF, {axisEta, axisPtHadron, axisMultFT0M, axisPosZ}, kAssocTrackNSteps});
    hAssocTracks->GetAxis(0)->SetTitle("#eta");
    hAssocTracks->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hAssocTracks->GetAxis(2)->SetTitle("multiplicity");
    hAssocTracks->GetAxis(3)->SetTitle("pos z");
  }

  /// Fill histograms of quantities independent from the daugther-mass hypothesis for data
  /// \param candidate is candidate
  template <typename T1>
  void fillHisto(const T1& candidate)
  {
    registry.fill(HIST("hPtCand"), candidate.pt());
    registry.fill(HIST("hEta"), candidate.eta());
    registry.fill(HIST("hEtaVsPtCand"), candidate.eta(), candidate.pt());
    registry.fill(HIST("hPhi"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));
    registry.fill(HIST("hPhiVsPtCand"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf), candidate.pt());
  }

  /// Fill histograms of quantities for the KKPi daugther-mass hypothesis for data
  /// \param candidate is candidate
  /// \param efficiencyWeight is the efficiency correction
  template <typename T1>
  void fillHistoKKPi(const T1& candidate, double efficiencyWeight)
  {
    registry.fill(HIST("hMassDsVsPt"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
    registry.fill(HIST("hMassDsData"), hfHelper.invMassDsToKKPi(candidate), efficiencyWeight);
    registry.fill(HIST("hSelectionStatusDsToKKPi"), candidate.isSelDsToKKPi());
  }

  /// Fill histograms of quantities for the PiKK daugther-mass hypothesis for data
  /// \param candidate is candidate
  /// \param efficiencyWeight is the efficiency correction
  template <typename T1>
  void fillHistoPiKK(const T1& candidate, double efficiencyWeight)
  {
    registry.fill(HIST("hMassDsVsPt"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
    registry.fill(HIST("hMassDsData"), hfHelper.invMassDsToPiKK(candidate), efficiencyWeight);
    registry.fill(HIST("hSelectionStatusDsToPiKK"), candidate.isSelDsToPiKK());
  }

  /// Fill histograms of quantities for the Ds signal for MC reco-level
  /// \param candidate is candidate
  /// \param multiplicityFT0M is the multiplicity
  template <typename T1>
  void fillHistoMcRecSig(const T1& candidate, float multiplicityFT0M)
  {
    registry.fill(HIST("hPtCandMcRecSig"), candidate.pt());
    registry.fill(HIST("hEtaMcRecSig"), candidate.eta());
    registry.fill(HIST("hPhiMcRecSig"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));

    // prompt and non-prompt division
    if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtCandMcRecSigPrompt"), candidate.pt());
      registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), candidate.pt(), multiplicityFT0M);
    } else if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
      registry.fill(HIST("hPtCandMcRecSigNonPrompt"), candidate.pt());
      registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), candidate.pt(), multiplicityFT0M);
    }
  }

  /// Fill histograms of quantities for the Ds backgroung for MC reco-level
  /// \param candidate is candidate
  template <typename T1>
  void fillHistoMcRecBkg(const T1& candidate)
  {
    registry.fill(HIST("hPtCandMcRecBkg"), candidate.pt());
    registry.fill(HIST("hEtaMcRecBkg"), candidate.eta());
    registry.fill(HIST("hPhiMcRecBkg"), RecoDecay::constrainAngle(candidate.phi(), -PIHalf));
  }

  /// Fill histograms of quantities for the Ds signal for MC reco-level
  /// \param particle is particle, Ds
  template <typename T1>
  void fillHistoMcGen(const T1& particle)
  {
    registry.fill(HIST("hPtCandMcGen"), particle.pt());
    registry.fill(HIST("hEtaMcGen"), particle.eta());
    registry.fill(HIST("hPhiMcGen"), RecoDecay::constrainAngle(particle.phi(), -PIHalf));

    // prompt and non-prompt division
    if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtCandMcGenPrompt"), particle.pt());
    } else if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
      registry.fill(HIST("hPtCandMcGenNonPrompt"), particle.pt());
    }
  }

  /// Ds-hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(SelCollisionsWithDs::iterator const& collision,
                   CandDsData const& candidates,
                   MyTracksData const& tracks)
  {
    BinningType corrBinning{{zPoolBins, multPoolBins}, true};
    registry.fill(HIST("hZVtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), collision.multFT0M());
    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    int nTracks = tracks.size();
    registry.fill(HIST("hMultiplicity"), nTracks);

    // Ds fill histograms and Ds-Hadron correlation for DsToKKPi
    for (const auto& candidate : candidates) {
      if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      double efficiencyWeightD = 1.;
      if (applyEfficiency) {
        efficiencyWeightD = 1. / efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt()));
      }
      std::vector<float> outputMl = {-1., -1., -1.};
      fillHisto(candidate);
      if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
        fillHistoKKPi(candidate, efficiencyWeightD);
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
        }
        entryDsCandRecoInfo(hfHelper.invMassDsToKKPi(candidate), candidate.pt(), outputMl[0], outputMl[2]);
      } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
        fillHistoPiKK(candidate, efficiencyWeightD);
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
        }
        entryDsCandRecoInfo(hfHelper.invMassDsToPiKK(candidate), candidate.pt(), outputMl[0], outputMl[2]);
      }
      if (candidate.isSelDsToKKPi() >= selectionFlagDs && candidate.isSelDsToPiKK() >= selectionFlagDs) {
        registry.fill(HIST("hCountSelectionStatusDsToKKPiAndToPiKK"), 0.);
      }

      // Ds-Hadron correlation dedicated section
      for (const auto& track : tracks) {
        // Removing Ds daughters by checking track indices
        if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
          continue;
        }
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }

        registry.fill(HIST("hEtaVsPtPartAssoc"), track.eta(), candidate.pt());
        registry.fill(HIST("hPhiVsPtPartAssoc"), RecoDecay::constrainAngle(track.phi(), -PIHalf), candidate.pt());
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(candidate), false);
          entryDsHadronGenInfo(false, false);
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(candidate), false);
          entryDsHadronGenInfo(false, false);
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
        }
      } // end track loop
    }   // end candidate loop
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processData, "Process data", true);

  /// Ds-Hadron correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(SelCollisionsWithDs::iterator const& collision,
                    CandDsMcReco const& candidates,
                    TracksWithMc const& tracks,
                    aod::McParticles const&)
  {
    BinningType corrBinning{{zPoolBins, multPoolBins}, true};
    registry.fill(HIST("hZVtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), collision.multFT0M());
    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    registry.fill(HIST("hCollisionPoolBin"), poolBin);

    // MC reco level
    bool isDsPrompt = false;
    bool isDsSignal = false;
    bool isAlreadyFilledEvent = false;
    float multiplicityFT0M = collision.multFT0M();
    for (const auto& candidate : candidates) {
      // prompt and non-prompt division
      isDsPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
      // Ds Signal
      isDsSignal = std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi;

      if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }

      double efficiencyWeightD = 1.;
      if (applyEfficiency) {
        efficiencyWeightD = 1. / efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt()));
      }
      if (isDsSignal) {
        fillHistoMcRecSig(candidate, multiplicityFT0M);
        // DsToKKPi and DsToPiKK division
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          registry.fill(HIST("hMassDsMcRec"), hfHelper.invMassDsToKKPi(candidate), efficiencyWeightD);
          registry.fill(HIST("hMassDsMcRecSig"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hMassDsVsPtMcRec"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hSelectionStatusDsToKKPi"), candidate.isSelDsToKKPi());
        }
        if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          registry.fill(HIST("hMassDsMcRec"), hfHelper.invMassDsToPiKK(candidate), efficiencyWeightD);
          registry.fill(HIST("hMassDsMcRecSig"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hMassDsVsPtMcRec"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hSelectionStatusDsToPiKK"), candidate.isSelDsToPiKK());
        }
        if (candidate.isSelDsToKKPi() >= selectionFlagDs && candidate.isSelDsToPiKK() >= selectionFlagDs) {
          registry.fill(HIST("hCountSelectionStatusDsToKKPiAndToPiKK"), 0.);
        }
      } else {
        fillHistoMcRecBkg(candidate);
        // DsToKKPi and DsToPiKK division
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          registry.fill(HIST("hMassDsMcRec"), hfHelper.invMassDsToKKPi(candidate), efficiencyWeightD);
          registry.fill(HIST("hMassDsMcRecBkg"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hMassDsVsPtMcRec"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hSelectionStatusDsToKKPi"), candidate.isSelDsToKKPi());
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          registry.fill(HIST("hMassDsMcRec"), hfHelper.invMassDsToPiKK(candidate), efficiencyWeightD);
          registry.fill(HIST("hMassDsMcRecBkg"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hMassDsVsPtMcRec"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeightD);
          registry.fill(HIST("hSelectionStatusDsToPiKK"), candidate.isSelDsToPiKK());
        }
      }
      std::vector<float> outputMl = {-1., -1., -1.};

      // Ds-Hadron correlation dedicated section
      // if the candidate is selected as Ds, search for Hadron and evaluate correlations
      for (const auto& track : tracks) {
        // Removing Ds daughters by checking track indices
        if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
          continue;
        }
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }
        bool isPhysicalPrimary = false;
        // DsToKKPi and DsToPiKK division
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(candidate), isDsSignal);
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
          }
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
          if (track.has_mcParticle()) {
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            isPhysicalPrimary = mcParticle.isPhysicalPrimary();
            entryDsHadronGenInfo(isDsPrompt, isPhysicalPrimary);
          } else {
            entryDsHadronGenInfo(isDsPrompt, isPhysicalPrimary);
            registry.fill(HIST("hFakeTracksMcRec"), track.pt());
          }
          // for secondary particle fraction estimation
          if (!isAlreadyFilledEvent) {
            registry.fill(HIST("hPtParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
            if (isPhysicalPrimary) {
              registry.fill(HIST("hPtPrimaryParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
            }
          }
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(candidate), isDsSignal);
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
          }
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
          if (track.has_mcParticle()) {
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            isPhysicalPrimary = mcParticle.isPhysicalPrimary();
            entryDsHadronGenInfo(isDsPrompt, isPhysicalPrimary);
          } else {
            entryDsHadronGenInfo(isDsPrompt, false);
            registry.fill(HIST("hFakeTracksMcRec"), track.pt());
          }
          // for secondary particle fraction estimation
          if (!isAlreadyFilledEvent) {
            registry.fill(HIST("hPtParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
            if (isPhysicalPrimary) {
              registry.fill(HIST("hPtPrimaryParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
            }
          }
        }
      } // end track loop
      isAlreadyFilledEvent = true;
    } // end candidate loop
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcRec, "Process MC Reco mode", false);

  /// Ds-Hadron correlation - for calculating candidate reconstruction efficiency using MC reco-level analysis
  void processMcCandEfficiency(soa::Join<aod::Collisions, aod::FT0Mults> const&,
                               soa::Join<aod::McCollisions, aod::MultsExtraMC> const&,
                               CandDsMcGen const& mcParticles,
                               CandDsMcReco const& candidates,
                               aod::TracksWMc const&)
  {
    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));

    /// Gen loop
    float multiplicity = -1.;
    for (auto& mcParticle : mcParticles) {
      // generated candidates
      if (std::abs(mcParticle.pdgCode()) == Pdg::kDS) {
        auto mcCollision = mcParticle.template mcCollision_as<soa::Join<aod::McCollisions, aod::MultsExtraMC>>();
        multiplicity = mcCollision.multMCFT0A() + mcCollision.multMCFT0C(); // multFT0M = multFt0A + multFT0C
        hCandidates->Fill(kCandidateStepMcGenAll, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
        if (std::abs(mcParticle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
          hCandidates->Fill(kCandidateStepMcGenDsToKKPi, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
          auto yDs = RecoDecay::y(mcParticle.pVector(), o2::constants::physics::MassDS);
          if (std::abs(yDs) <= yCandGenMax) {
            hCandidates->Fill(kCandidateStepMcCandInAcceptance, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
          }
          bool isDaughterInAcceptance = true;
          auto daughters = mcParticle.template daughters_as<CandDsMcGen>();
          for (const auto& daughter : daughters) {
            if (daughter.pt() < ptDaughterMin || std::abs(daughter.eta()) > etaTrackMax) {
              isDaughterInAcceptance = false;
            }
          }
          if (isDaughterInAcceptance) {
            hCandidates->Fill(kCandidateStepMcDaughtersInAcceptance, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
            fillHistoMcGen(mcParticle);
          }
        }
      }
    }

    // recontructed candidates loop
    for (auto& candidate : candidates) {
      auto collision = candidate.template collision_as<soa::Join<aod::Collisions, aod::FT0Mults>>();
      multiplicity = collision.multFT0M();
      registry.fill(HIST("hMultFT0M"), multiplicity);
      if (std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        auto prong0McPart = candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<CandDsMcGen>();
        // DsToKKPi and DsToPiKK division
        if (((std::abs(prong0McPart.pdgCode()) == kKPlus) && (candidate.isSelDsToKKPi() >= selectionFlagDs)) || ((std::abs(prong0McPart.pdgCode()) == kPiPlus) && (candidate.isSelDsToPiKK() >= selectionFlagDs))) {
          registry.fill(HIST("hPtCand"), candidate.pt());
          hCandidates->Fill(kCandidateStepMcReco, candidate.pt(), multiplicity, candidate.originMcRec());
          if (std::abs(hfHelper.yDs(candidate)) <= yCandMax) {
            hCandidates->Fill(kCandidateStepMcRecoInAcceptance, candidate.pt(), multiplicity, candidate.originMcRec());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcCandEfficiency, "Process MC for calculating candidate reconstruction efficiency", false);

  /// Ds-Hadron correlation - for calculating associated particle tracking efficiency using MC reco-level analysis
  void processMcTrackEfficiency(soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels> const&,
                                soa::Join<aod::McCollisions, aod::MultsExtraMC> const&,
                                CandDsMcGen const& mcParticles,
                                TracksWithMc const& tracksData)
  {
    auto hAssocTracks = registry.get<StepTHn>(HIST("hAssocTracks"));

    /// Gen loop
    float multiplicity = -1.;
    float posZ = -20.;
    for (auto& mcParticle : mcParticles) {
      // generated tracks
      if (mcParticle.isPhysicalPrimary() && ((std::abs(mcParticle.pdgCode()) == kElectron) || (std::abs(mcParticle.pdgCode()) == kMuonMinus) || (std::abs(mcParticle.pdgCode()) == kPiPlus) || (std::abs(mcParticle.pdgCode()) == kKPlus) || (std::abs(mcParticle.pdgCode()) == kProton))) {
        auto mcCollision = mcParticle.template mcCollision_as<soa::Join<aod::McCollisions, aod::MultsExtraMC>>();
        multiplicity = mcCollision.multMCFT0A() + mcCollision.multMCFT0C(); // multFT0M = multFt0A + multFT0C
        posZ = mcCollision.posZ();
        hAssocTracks->Fill(kAssocTrackStepMcGen, mcParticle.eta(), mcParticle.pt(), multiplicity, posZ);
        if (mcParticle.pt() > ptTrackMin && std::abs(mcParticle.eta()) < etaTrackMax) {
          hAssocTracks->Fill(kAssocTrackStepMcGenInAcceptance, mcParticle.eta(), mcParticle.pt(), multiplicity, posZ);
          registry.fill(HIST("hPtParticleAssocMcGen"), mcParticle.pt());
        }
      }
    }

    // recontructed tracks loop
    for (auto& track : tracksData) {
      if (track.has_collision()) {
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }
        auto collision = track.template collision_as<soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels>>();
        if (useSel8ForTrackEff && !collision.sel8()) {
          continue;
        }
        multiplicity = collision.multFT0M();
        posZ = collision.posZ();
        registry.fill(HIST("hZVtx"), posZ);
        registry.fill(HIST("hMultFT0M"), multiplicity);
        hAssocTracks->Fill(kAssocTrackStepRecoAll, track.eta(), track.pt(), multiplicity, posZ);
        if (track.has_mcParticle()) {
          auto mcParticle = track.template mcParticle_as<CandDsMcGen>();
          hAssocTracks->Fill(kAssocTrackStepRecoMcMatch, mcParticle.eta(), mcParticle.pt(), multiplicity, posZ);
          if (mcParticle.isPhysicalPrimary()) {
            hAssocTracks->Fill(kAssocTrackStepRecoPrimaries, mcParticle.eta(), mcParticle.pt(), multiplicity, posZ);
            registry.fill(HIST("hPtParticleAssocMcRec"), track.pt());
          }
        }
      } else {
        // fake track
        hAssocTracks->Fill(kAssocTrackStepFake, track.eta(), track.pt(), multiplicity, posZ);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcTrackEfficiency, "Process MC for calculating associated particle tracking efficiency", false);

  /// Ds-Hadron correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(SelCollisionsWithDsMc::iterator const& mcCollision,
                    CandDsMcGen const& mcParticles)
  {
    auto getTracksSize = [&mcParticles](SelCollisionsWithDsMc::iterator const&) {
      int nTracks = 0;
      for (const auto& track : mcParticles) {
        if (track.isPhysicalPrimary() && std::abs(track.eta()) < 1.0) {
          nTracks++;
        }
      }
      return nTracks;
    };
    using BinningTypeMcGen = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::mccollision::PosZ, decltype(getTracksSize)>;
    BinningTypeMcGen corrBinningMcGen{{getTracksSize}, {zPoolBins, multPoolBins}, true};                    // TODO
    int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), getTracksSize(mcCollision))); // TODO
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    bool isDsPrompt = false;

    // MC gen level
    for (const auto& particle : mcParticles) {
      // check if the particle is Ds
      if (std::abs(particle.pdgCode()) != Pdg::kDS) {
        continue;
      }
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        double yD = RecoDecay::y(particle.pVector(), MassDS); // TODO
        if (yCandGenMax >= 0. && std::abs(yD) > yCandGenMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
          continue;
        }
        fillHistoMcGen(particle);
        // prompt and non-prompt division
        isDsPrompt = particle.originMcGen() == RecoDecay::OriginType::Prompt;

        // Ds Hadron correlation dedicated section
        for (const auto& particleAssoc : mcParticles) {
          if (std::abs(particleAssoc.eta()) > etaTrackMax) {
            continue;
          }
          if (particleAssoc.pt() < ptTrackMin) {
            continue;
          }
          if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particleAssoc.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
            continue;
          }
          registry.fill(HIST("hPtParticleAssocMcGen"), particleAssoc.pt());
          entryDsHadronPair(getDeltaPhi(particleAssoc.phi(), particle.phi()),
                            particleAssoc.eta() - particle.eta(),
                            particle.pt(),
                            particleAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(MassDS, true);
          entryDsHadronGenInfo(isDsPrompt, particleAssoc.isPhysicalPrimary());
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcGen, "Process MC Gen mode", false);

  // Event Mixing
  void processDataME(SelCollisionsWithDs const& collisions,
                     CandDsData const& candidates,
                     MyTracksData const& tracks)
  {
    BinningType corrBinning{{zPoolBins, multPoolBins}, true};
    for (const auto& collision : collisions) {
      registry.fill(HIST("hMultFT0M"), collision.multFT0M());
      registry.fill(HIST("hZVtx"), collision.posZ());
    }

    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelCollisionsWithDs, CandDsData, MyTracksData, BinningType> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      if (tracks1.size() == 0) {
        continue;
      }

      // LOGF(info, "Mixed event collisions: Index = (%d, %d), tracks Size: (%d, %d), Z Vertex: (%f, %f), Pool Bin: (%d, %d)", c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size(), c1.posZ(), c2.posZ(), corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M())), corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M())));
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      int poolBinDs = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M()));
      registry.fill(HIST("hTracksPoolBin"), poolBin); // note that the selections here are not yet applied
      registry.fill(HIST("hDsPoolBin"), poolBinDs);   // note that the selections here are not yet applied
      for (const auto& [cand, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!(cand.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) {
          continue;
        }
        if (std::abs(hfHelper.yDs(cand)) > yCandMax || cand.pt() < ptCandMin || cand.pt() > ptCandMax) {
          continue;
        }
        if (!pAssoc.isGlobalTrackWoDCA()) {
          continue;
        }
        std::vector<float> outputMl = {-1., -1., -1.};
        // DsToKKPi and DsToPiKK division
        if (cand.isSelDsToKKPi() >= selectionFlagDs) {
          // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d), KKPi", cand.index(), pAssoc.index(), c1.index(), c2.index(), cand.collision().index(), pAssoc.collision().index());
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), cand.phi()),
                            pAssoc.eta() - cand.eta(),
                            cand.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(cand), false);
          entryDsHadronGenInfo(false, false);
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = cand.mlProbDsToKKPi()[classMl->at(iclass)];
          }
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
        } else if (cand.isSelDsToPiKK() >= selectionFlagDs) {
          // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d), PiKK", cand.index(), pAssoc.index(), c1.index(), c2.index(), cand.collision().index(), pAssoc.collision().index());
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), cand.phi()),
                            pAssoc.eta() - cand.eta(),
                            cand.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(cand), false);
          entryDsHadronGenInfo(false, false);
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = cand.mlProbDsToPiKK()[classMl->at(iclass)];
          }
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processDataME, "Process Mixed Event Data", false);

  void processMcRecME(SelCollisionsWithDs const& collisions,
                      CandDsMcReco const& candidates,
                      MyTracksData const& tracks)
  {
    BinningType corrBinning{{zPoolBins, multPoolBins}, true};
    for (const auto& candidate : candidates) {
      if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      if (std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        // DsToKKPi and DsToPiKK division
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          fillHistoMcRecSig(candidate, 0.);
          registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelDsToKKPi());
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          fillHistoMcRecSig(candidate, 0.);
          registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelDsToPiKK());
        }
      } else {
        fillHistoMcRecBkg(candidate);
      }
    }
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelCollisionsWithDs, CandDsMcReco, MyTracksData, BinningType> pairMcRec{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    bool isDsPrompt = false;
    bool isDsSignal = false;
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcRec) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      int poolBinDs = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M()));
      registry.fill(HIST("hMultFT0M"), c1.multFT0M());
      registry.fill(HIST("hZVtx"), c1.posZ());
      registry.fill(HIST("hTracksPoolBin"), poolBin); // note that the selections here are not yet applied
      registry.fill(HIST("hDsPoolBin"), poolBinDs);   // note that the selections here are not yet applied
      for (const auto& [candidate, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (std::abs(hfHelper.yDs(candidate)) > yCandMax || candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
          continue;
        }
        if (!pAssoc.isGlobalTrackWoDCA()) {
          continue;
        }
        std::vector<float> outputMl = {-1., -1., -1.};
        // prompt and non-prompt division
        isDsPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
        // Ds Signal
        isDsSignal = std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi;
        // DsToKKPi and DsToPiKK division
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), candidate.phi()),
                            pAssoc.eta() - candidate.eta(),
                            candidate.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(candidate), isDsSignal);
          entryDsHadronGenInfo(isDsPrompt, false);
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
          }
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), candidate.phi()),
                            pAssoc.eta() - candidate.eta(),
                            candidate.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(candidate), isDsSignal);
          entryDsHadronGenInfo(isDsPrompt, false);
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
          }
          entryDsHadronMlInfo(outputMl[0], outputMl[2]);
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcRecME, "Process Mixed Event MCRec", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDsHadronsSelCollision>(cfgc),
                      adaptAnalysisTask<HfCorrelatorDsHadrons>(cfgc)};
}
