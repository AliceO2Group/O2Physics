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
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
double getDeltaPhi(double phiD, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PIHalf);
}

/// definition of variables for Ds hadron pairs (in data-like, MC-reco and MC-kine tasks)
const int nBinsPtMassAndEfficiency = o2::analysis::hf_cuts_ds_to_k_k_pi::nBinsPt;
const double efficiencyDmesonDefault[nBinsPtMassAndEfficiency] = {};
auto vecEfficiencyDmeson = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + nBinsPtMassAndEfficiency};

// histograms axes definition
AxisSpec axisMassD = {200, 1.75, 2.13, "inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
AxisSpec axisPhi = {128, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf, "#it{#varphi}"};
AxisSpec axisEta = {100, -2., 2., "#it{#eta}"};
AxisSpec axisY = {100, -2., 2., "#it{#y}"};
AxisSpec axisPtD = {180, 0., 36., "#it{p}_{T}Ds (GeV/#it{c})"};
AxisSpec axisPtProng0 = {180, 0., 36., "#it{p}_{T} prong0 (GeV/#it{c})"};
AxisSpec axisPtProng1 = {180, 0., 36., "#it{p}_{T} prong1 (GeV/#it{c})"};
AxisSpec axisPtProng2 = {180, 0., 36., "#it{p}_{T} prong2 (GeV/#it{c})"};
AxisSpec axisPtHadron = {11, 0., 11., "#it{p}_{T} Hadron (GeV/#it{c})"};
AxisSpec axisMultiplicity = {1000, 0., 10000., "Multiplicity"};
AxisSpec axisPoolBin = {9, 0., 9., "PoolBin"};

// binning type
std::vector<double> zBins{VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0};
// std::vector<double> zBins{VARIABLE_WIDTH, -10.0, 10.0};
std::vector<double> multBins{VARIABLE_WIDTH, 0., 500.0, 5000., 10000.};
// std::vector<double> multBins{VARIABLE_WIDTH, 0., 500.0, 5000.};
std::vector<double> multBinsMcGen{VARIABLE_WIDTH, 0., 20., 50.0, 500.}; // In MCGen multiplicity is defined by counting primaries
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
BinningType corrBinning{{zBins, multBins}, true};

/// Code to select collisions with at least one Ds meson
struct HfCorrelatorDsHadronsSelCollision {
  Produces<aod::DmesonSelection> collisionsWithSelDs;

  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  SliceCache cache;
  HfHelper hfHelper;

  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Filter dsFilter = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) != static_cast<uint8_t>(0); // filter in HfCand3Prong

  Partition<CandDsData> selectedDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;
  Partition<CandDsMcReco> selectedDsCandidatesMc = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;

  /// Code to select collisions with at least one Ds meson - for real data and data-like analysis
  void processDsSelCollisionsData(aod::Collision const& collision,
                                  CandDsData const& candidates)
  {
    bool isDsFound = false;
    if (selectedDsCandidates.size() > 0) {
      auto selectedDsAllCandGrouped = selectedDsCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      for (const auto& candidate : selectedDsAllCandGrouped) {
        if (yCandMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
          continue;
        }
        isDsFound = true;
        break;
      }
    }
    collisionsWithSelDs(isDsFound);
  }
  PROCESS_SWITCH(HfCorrelatorDsHadronsSelCollision, processDsSelCollisionsData, "Process Ds Collision Selection Data", true);

  /// Code to select collisions with at least one Ds meson - for MC reco-level analysis
  void processDsSelCollisionsMcRec(aod::Collision const& collision,
                                   CandDsMcReco const& candidates)
  {
    bool isDsFound = false;
    if (selectedDsCandidatesMc.size() > 0) {
      auto selectedDsCandidatesGroupedMc = selectedDsCandidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      for (const auto& candidate : selectedDsCandidatesGroupedMc) {
        if (yCandMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
          continue;
        }
        isDsFound = true;
        break;
      }
    }
    collisionsWithSelDs(isDsFound);
  }
  PROCESS_SWITCH(HfCorrelatorDsHadronsSelCollision, processDsSelCollisionsMcRec, "Process Ds Collision Selection MCRec", false);

  /// Code to select collisions with at least one Ds meson - for MC gen-level analysis
  void processDsSelCollisionsMcGen(aod::McCollision const& mcCollision,
                                   CandDsMcGen const& mcParticles)
  {
    bool isDsFound = false;
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != pdg::Code::kDS) {
        continue;
      }
      double yD = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::analysis::pdg::MassDS);
      if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
        continue;
      }
      isDsFound = true;
      break;
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

  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<bool> applyEfficiency{"applyEfficiency", true, "Flag for applying D-meson efficiency weights"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCA_xy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCA_z of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "max. cand pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 50., "max. track pT"};
  Configurable<float> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<float> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{vecEfficiencyDmeson}, "Efficiency values for Ds meson"};

  HfHelper hfHelper;
  SliceCache cache;

  using SelCollisionsWithDs = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::DmesonSelection>>;      // collisionFilter applied
  using SelCollisionsWithDsMc = soa::Filtered<soa::Join<aod::McCollisions, aod::DmesonSelection>>;              // collisionFilter applied
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;                           // flagDsFilter applied
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>; // flagDsFilter applied
  using CandDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;                                      // flagDsFilter applied
  using MyTracksData = soa::Filtered<aod::TracksWDca>;                                                          // trackFilter applied

  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter flagDsFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs);
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);

  Preslice<aod::HfCand3Prong> perCol = aod::hf_cand::collisionId;

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "Ds,Hadron candidates", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng0", "Ds,Hadron candidates - prong 0", {HistType::kTH1F, {axisPtProng0}}},
     {"hPtProng1", "Ds,Hadron candidates - prong 1", {HistType::kTH1F, {axisPtProng1}}},
     {"hPtProng2", "Ds,Hadron candidates - prong 2", {HistType::kTH1F, {axisPtProng2}}},
     {"hSelectionStatusDsToKKPi", "Ds,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hSelectionStatusDsToPiKK", "Ds,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hEta", "Ds,Hadron candidates", {HistType::kTH1F, {axisEta}}},
     {"hEtaVsPtCand", "Ds,Hadron candidates", {HistType::kTH2F, {{axisEta}, {axisPtD}}}},
     {"hEtaVsPtPartAssoc", "Particles associated", {HistType::kTH2F, {{axisEta}, {axisPtD}}}},
     {"hPhi", "Ds,Hadron candidates", {HistType::kTH1F, {axisPhi}}},
     {"hPhiVsPtCand", "Ds,Hadron candidates", {HistType::kTH2F, {{axisPhi}, {axisPtD}}}},
     {"hPhiVsPtPartAssoc", "Particles associated", {HistType::kTH2F, {{axisPhi}, {axisPtD}}}},
     {"hY", "Ds,Hadron candidates", {HistType::kTH1F, {axisY}}},
     {"hPtCandMcRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtD}}},
     {"hPtCandMcRecSigPrompt", "Ds,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtD}}},
     {"hPtCandMcRecSigNonPrompt", "Ds,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng0McRecSig", "Ds,Hadron candidates - MC Reco - prong 0", {HistType::kTH1F, {axisPtProng0}}},
     {"hPtProng1McRecSig", "Ds,Hadron candidates - MC Reco - prong 1", {HistType::kTH1F, {axisPtProng1}}},
     {"hPtProng2McRecSig", "Ds,Hadron candidates - MC Reco - prong 2", {HistType::kTH1F, {axisPtProng2}}},
     {"hPtCandMcRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng0McRecBkg", "Ds,Hadron candidates - MC Reco - prong 0", {HistType::kTH1F, {axisPtProng0}}},
     {"hPtProng1McRecBkg", "Ds,Hadron candidates - MC Reco - prong 1", {HistType::kTH1F, {axisPtProng1}}},
     {"hPtProng2McRecBkg", "Ds,Hadron candidates - MC Reco - prong 2", {HistType::kTH1F, {axisPtProng2}}},
     {"hSelectionStatusMcRec", "Ds,Hadron candidates - MC Reco;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hEtaMcRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}}},
     {"hPhiMcRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}}},
     {"hYMcRecSig", "Ds,Hadron candidates - MC Reco;;entries", {HistType::kTH1F, {axisY}}},
     {"hEtaMcRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}}},
     {"hPhiMcRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}}},
     {"hYMcRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisY}}},
     {"hMcEvtCount", "Event counter - MC Gen", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMcGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisPtD}}},
     {"hPtCandMcGenPrompt", "Ds,Hadron particles - MC Gen Prompt", {HistType::kTH1F, {axisPtD}}},
     {"hPtCandMcGenNonPrompt", "Ds,Hadron particles - MC Gen Non Prompt", {HistType::kTH1F, {axisPtD}}},
     {"hPtParticleAssocMcRec", "Associated Particle - MC Rec", {HistType::kTH1F, {axisPtHadron}}},
     {"hPtParticleAssocMcGen", "Associated Particle - MC Gen", {HistType::kTH1F, {axisPtHadron}}},
     {"hEtaMcGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisEta}}},
     {"hPhiMcGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisPhi}}},
     {"hYMcGen", "Ds,Hadron candidates - MC Gen", {HistType::kTH1F, {axisY}}},
     {"hCountDsHadronPerEvent", "Ds,Hadron particles - MC Gen;Number per event;entries", {HistType::kTH1F, {{21, -0.5, 20.5}}}},
     {"hMultiplicityPreSelection", "Multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {axisMultiplicity}}},
     {"hPtVsMultiplicityMcRecPrompt", "Multiplicity V0M - MC Rec Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultiplicity}}}},
     {"hPtVsMultiplicityMcRecNonPrompt", "Multiplicity V0M - MC Rec Non Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultiplicity}}}},
     {"hPtVsMultiplicityMcGenPrompt", "Multiplicity V0M - MC Gen Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultiplicity}}}},
     {"hPtVsMultiplicityMcGenNonPrompt", "Multiplicity V0M - MC Gen Non Prompt", {HistType::kTH2F, {{axisPtD}, {axisMultiplicity}}}},
     {"hMultiplicity", "Multiplicity", {HistType::kTH1F, {axisMultiplicity}}},
     {"hMultV0M", "multiplicity;multiplicity;entries", {HistType::kTH1F, {axisMultiplicity}}},
     {"hZVtx", "z vertex;z vertex;entries", {HistType::kTH1F, {{200, -20., 20.}}}},
     {"hCollisionPoolBin", "Collisions in each pool Bin;pool Bin;entries", {HistType::kTH1F, {axisPoolBin}}},
     {"hDsPoolBin", "Ds selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {axisPoolBin}}},
     {"hTracksPoolBin", "Tracks selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {axisPoolBin}}}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMassDsVsPt", "Ds candidates", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDsData", "Ds candidates", {HistType::kTH1F, {axisMassD}});
    registry.add("hMassDsMCRec", "Ds candidates", {HistType::kTH1F, {axisMassD}});
    registry.add("hMassDsVsPtMCRec", "Ds signal candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDsMCRecSig", "Ds signal candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDsMCRecBkg", "Ds background candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCountDstriggersMCGen", "Ds trigger particles - MC Gen", {HistType::kTH2F, {{1, -0.5, 0.5, "number of Ds triggers"}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  /// Fill histograms of quantities independent from the daugther-mass hypothesis for data
  /// \param candidate is candidate
  template <typename T1>
  void fillHisto(const T1& candidate)
  {
    registry.fill(HIST("hPtCand"), candidate.pt());
    registry.fill(HIST("hPtProng0"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1"), candidate.ptProng1());
    registry.fill(HIST("hPtProng2"), candidate.ptProng2());
    registry.fill(HIST("hEta"), candidate.eta());
    registry.fill(HIST("hEtaVsPtCand"), candidate.eta(), candidate.pt());
    registry.fill(HIST("hPhi"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf));
    registry.fill(HIST("hPhiVsPtCand"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf), candidate.pt());
    registry.fill(HIST("hY"), hfHelper.yDs(candidate));
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
  /// \param multiplicityV0M is the multiplicity
  template <typename T1>
  void fillHistoMcRecSig(const T1& candidate, float multiplicityV0M)
  {
    registry.fill(HIST("hPtCandMcRecSig"), candidate.pt());
    registry.fill(HIST("hPtProng0McRecSig"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1McRecSig"), candidate.ptProng1());
    registry.fill(HIST("hPtProng2McRecSig"), candidate.ptProng2());
    registry.fill(HIST("hEtaMcRecSig"), candidate.eta());
    registry.fill(HIST("hPhiMcRecSig"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf));
    registry.fill(HIST("hYMcRecSig"), hfHelper.yDs(candidate));

    // prompt and non-prompt division
    if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtCandMcRecSigPrompt"), candidate.pt());
      registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), candidate.pt(), multiplicityV0M);
    } else if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
      registry.fill(HIST("hPtCandMcRecSigNonPrompt"), candidate.pt());
      registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), candidate.pt(), multiplicityV0M);
    }
  }

  /// Fill histograms of quantities for the Ds backgroung for MC reco-level
  /// \param candidate is candidate
  template <typename T1>
  void fillHistoMcRecBkg(const T1& candidate)
  {
    registry.fill(HIST("hPtCandMcRecBkg"), candidate.pt());
    registry.fill(HIST("hPtProng0McRecBkg"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1McRecBkg"), candidate.ptProng1());
    registry.fill(HIST("hPtProng2McRecBkg"), candidate.ptProng2());
    registry.fill(HIST("hEtaMcRecBkg"), candidate.eta());
    registry.fill(HIST("hPhiMcRecBkg"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf));
    registry.fill(HIST("hYMcRecBkg"), hfHelper.yDs(candidate));
  }

  /// Fill histograms of quantities for the Ds signal for MC reco-level
  /// \param particle is particle, Ds
  /// \param yD is the Ds rapidity
  template <typename T1>
  void fillHistoMcGen(const T1& particle, double yD)
  {
    registry.fill(HIST("hPtCandMcGen"), particle.pt());
    registry.fill(HIST("hEtaMcGen"), particle.eta());
    registry.fill(HIST("hPhiMcGen"), RecoDecay::constrainAngle(particle.phi(), -o2::constants::math::PIHalf));
    registry.fill(HIST("hYMcGen"), yD);
    registry.fill(HIST("hCountDstriggersMCGen"), 0, particle.pt()); // to count trigger Ds for normalisation

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
    // if (selectedDsCandidates.size() > 0) {
    if (candidates.size() > 0) {
      registry.fill(HIST("hZVtx"), collision.posZ());
      registry.fill(HIST("hMultV0M"), collision.multFV0M());
      int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFV0M()));
      registry.fill(HIST("hCollisionPoolBin"), poolBin);
      int nTracks = 0;
      if (collision.numContrib() > 1) {
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          nTracks++;
          registry.fill(HIST("hTracksPoolBin"), poolBin);
        }
      }
      if (nTracks < multMin || nTracks > multMax) {
        return;
      }
      registry.fill(HIST("hMultiplicity"), nTracks);

      // Ds fill histograms and Ds-Hadron correlation for DsToKKPi
      for (const auto& candidate : candidates) {
        if (yCandMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
          continue;
        }
        if (candidate.pt() > ptCandMax) {
          continue;
        }
        double efficiencyWeight = 1.;
        if (applyEfficiency) {
          efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate.pt()));
        }
        fillHisto(candidate);
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          fillHistoKKPi(candidate, efficiencyWeight);
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          fillHistoPiKK(candidate, efficiencyWeight);
        }

        // Ds-Hadron correlation dedicated section
        for (const auto& track : tracks) {
          // Removing Ds daughters by checking track indices
          if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
            continue;
          }
          registry.fill(HIST("hEtaVsPtPartAssoc"), track.eta(), candidate.pt());
          registry.fill(HIST("hPhiVsPtPartAssoc"), RecoDecay::constrainAngle(track.phi(), -o2::constants::math::PIHalf), candidate.pt());
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
            entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                              track.eta() - candidate.eta(),
                              candidate.pt(),
                              track.pt(),
                              poolBin);
            entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(candidate), false);
            entryDsHadronGenInfo(false);
          } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
            entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                              track.eta() - candidate.eta(),
                              candidate.pt(),
                              track.pt(),
                              poolBin);
            entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(candidate), false);
            entryDsHadronGenInfo(false);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processData, "Process data", true);

  /// Ds-Hadron correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(SelCollisionsWithDs::iterator const& collision,
                    CandDsMcReco const& candidates,
                    MyTracksData const& tracks)
  {
    if (candidates.size() > 0) {
      registry.fill(HIST("hZVtx"), collision.posZ());
      registry.fill(HIST("hMultV0M"), collision.multFV0M());
      int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFV0M()));
      registry.fill(HIST("hCollisionPoolBin"), poolBin);
      int nTracks = 0;
      if (collision.numContrib() > 1) {
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          registry.fill(HIST("hTracksPoolBin"), poolBin);
          nTracks++;
        }
      }
      registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
      if (nTracks < multMin || nTracks > multMax) {
        return;
      }
      registry.fill(HIST("hMultiplicity"), nTracks);

      // auto selectedDsMcRecoCandGrouped = selectedDsMcRecoCand->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

      // MC reco level
      bool isDsPrompt = false;
      bool isDsSignal = false;
      float multiplicityV0M = collision.multFV0M();
      for (const auto& candidate : candidates) {
        // prompt and non-prompt division
        isDsPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
        // Ds Signal
        isDsSignal = std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi;

        if (yCandMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
          continue;
        }
        if (candidate.pt() >= ptCandMax) {
          continue;
        }
        double efficiencyWeight = 1.;
        if (applyEfficiency) {
          efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate.pt()));
        }
        if (isDsSignal) {
          fillHistoMcRecSig(candidate, multiplicityV0M);
          // DsToKKPi and DsToPiKK division
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
            registry.fill(HIST("hMassDsMCRec"), hfHelper.invMassDsToKKPi(candidate), efficiencyWeight);
            registry.fill(HIST("hMassDsMCRecSig"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hMassDsVsPtMCRec"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelDsToKKPi());
          } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
            registry.fill(HIST("hMassDsMCRec"), hfHelper.invMassDsToPiKK(candidate), efficiencyWeight);
            registry.fill(HIST("hMassDsMCRecSig"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hMassDsVsPtMCRec"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelDsToPiKK());
          }
        } else {
          fillHistoMcRecBkg(candidate);
          // DsToKKPi and DsToPiKK division
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
            registry.fill(HIST("hMassDsMCRec"), hfHelper.invMassDsToKKPi(candidate), efficiencyWeight);
            registry.fill(HIST("hMassDsMCRecBkg"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hMassDsVsPtMCRec"), hfHelper.invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelDsToKKPi());
          } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
            registry.fill(HIST("hMassDsMCRec"), hfHelper.invMassDsToPiKK(candidate), efficiencyWeight);
            registry.fill(HIST("hMassDsMCRecBkg"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hMassDsVsPtMCRec"), hfHelper.invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelDsToPiKK());
          }
        }

        // Ds-Hadron correlation dedicated section
        // if the candidate is selected as Ds, search for Hadron and evaluate correlations
        for (const auto& track : tracks) {
          // Removing Ds daughters by checking track indices
          if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
            continue;
          }
          registry.fill(HIST("hPtParticleAssocMcRec"), track.pt()); // va tolto
          // DsToKKPi and DsToPiKK division
          if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
            entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                              track.eta() - candidate.eta(),
                              candidate.pt(),
                              track.pt(),
                              poolBin);
            entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(candidate), isDsSignal);
            entryDsHadronGenInfo(isDsPrompt);
          } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
            entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                              track.eta() - candidate.eta(),
                              candidate.pt(),
                              track.pt(),
                              poolBin);
            entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(candidate), isDsSignal);
            entryDsHadronGenInfo(isDsPrompt);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcRec, "Process MC Reco mode", false);

  /// Ds-Hadron correlation - for calculating efficiencies using MC reco-level analysis
  void processMcEfficiencies(CandDsMcReco const& candidates,
                             CandDsMcGen const& mcParticles,
                             MyTracksData const& tracksData,
                             aod::TracksWMc const&)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      if (yCandMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
        continue;
      }
      if (candidate.pt() > ptCandMax) {
        continue;
      }
      if (std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        auto prong0McPart = candidate.prong0_as<aod::TracksWMc>().mcParticle_as<CandDsMcGen>();
        // DsToKKPi and DsToPiKK division
        if ((std::abs(prong0McPart.pdgCode()) == kKPlus) && (candidate.isSelDsToKKPi() >= selectionFlagDs)) {
          fillHistoMcRecSig(candidate, 0.);
          registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelDsToKKPi());
        } else if ((std::abs(prong0McPart.pdgCode()) == kPiPlus) && (candidate.isSelDsToPiKK() >= selectionFlagDs)) {
          fillHistoMcRecSig(candidate, 0.);
          registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelDsToPiKK());
        }
      } else {
        fillHistoMcRecBkg(candidate);
      }
    }

    // MC reco level for particles associated reconstruction's efficiency
    for (const auto& track : tracksData) {
      if (std::abs(track.eta()) > etaTrackMax) {
        continue;
      }
      if (track.pt() < ptTrackMin) {
        continue;
      }
      // Remove secondary tracks
      if (std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
        continue;
      }
      registry.fill(HIST("hPtParticleAssocMcRec"), track.pt());
    }

    // MC gen level for Ds meson reconstruction's efficiency
    for (const auto& particle : mcParticles) {
      // check if the particle is Ds
      if (std::abs(particle.pdgCode()) != pdg::Code::kDS) {
        continue;
      }
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        double yD = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::analysis::pdg::MassDS);
        if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
          continue;
        }
        fillHistoMcGen(particle, yD);
      }
    }

    // MC gen level for particles associated reconstruction's efficiency
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
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcEfficiencies, "Process MC for calculating efficiencies", false);

  /// Ds-Hadron correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(SelCollisionsWithDsMc::iterator const& mcCollision,
                    CandDsMcGen const& mcParticles)
  {
    int counterDsHadron = 0;
    registry.fill(HIST("hMcEvtCount"), 0);

    auto getTracksSize = [&mcParticles](SelCollisionsWithDsMc::iterator const& mcCollision) {
      int nTracks = 0;
      for (const auto& track : mcParticles) {
        if (track.isPhysicalPrimary() && std::abs(track.eta()) < 1.0) {
          nTracks++;
        }
      }
      return nTracks;
    };
    using BinningTypeMcGen = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::mccollision::PosZ, decltype(getTracksSize)>;
    BinningTypeMcGen corrBinningMcGen{{getTracksSize}, {zBins, multBinsMcGen}, true};
    int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), getTracksSize(mcCollision)));
    registry.fill(HIST("hCollisionPoolBin"), poolBin);
    bool isDsPrompt = false;

    // MC gen level
    for (const auto& particle : mcParticles) {
      // check if the particle is Ds
      if (std::abs(particle.pdgCode()) != pdg::Code::kDS) {
        continue;
      }
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DsToKKPi) {
        double yD = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::analysis::pdg::MassDS);
        if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
          continue;
        }
        fillHistoMcGen(particle, yD);
        counterDsHadron++;
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
          entryDsHadronRecoInfo(o2::analysis::pdg::MassDS, true);
          entryDsHadronGenInfo(isDsPrompt);
        }
      }
    }
    registry.fill(HIST("hCountDsHadronPerEvent"), counterDsHadron);
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcGen, "Process MC Gen mode", false);

  // Event Mixing
  void processDataME(soa::Join<aod::Collisions, aod::Mults> const& collisions,
                     CandDsData const& candidates,
                     MyTracksData const& tracks)
  {
    if (candidates.size() == 0) {
      return;
    }
    for (const auto& candidate : candidates) {
      fillHisto(candidate);
    }

    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<soa::Join<aod::Collisions, aod::Mults>, CandDsData, MyTracksData, BinningType> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      if (tracks1.size() == 0) {
        continue;
      }
      // LOGF(info, "Mixed event collisions: Index = (%d, %d), tracks Size: (%d, %d), Z Vertex: (%f, %f), Pool Bin: (%d, %d)", c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size(), c1.posZ(), c2.posZ(), corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFV0M())), corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFV0M())));
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFV0M()));
      int poolBinDs = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFV0M()));
      registry.fill(HIST("hTracksPoolBin"), poolBin);
      registry.fill(HIST("hDsPoolBin"), poolBinDs);
      for (const auto& [cand, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!(cand.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(hfHelper.yDs(cand)) > yCandMax) {
          continue;
        }

        // DsToKKPi and DsToPiKK division
        if (cand.isSelDsToKKPi() >= selectionFlagDs) {
          // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d), KKPi", cand.index(), pAssoc.index(), c1.index(), c2.index(), cand.collision().index(), pAssoc.collision().index());
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), cand.phi()),
                            pAssoc.eta() - cand.eta(),
                            cand.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToKKPi(cand), false);
          entryDsHadronGenInfo(false);
        } else if (cand.isSelDsToPiKK() >= selectionFlagDs) {
          // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d), PiKK", cand.index(), pAssoc.index(), c1.index(), c2.index(), cand.collision().index(), pAssoc.collision().index());
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), cand.phi()),
                            pAssoc.eta() - cand.eta(),
                            cand.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(cand), false);
          entryDsHadronGenInfo(false);
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processDataME, "Process Mixed Event Data", false);

  void processMcRecME(SelCollisionsWithDs const& collisions,
                      CandDsMcReco const& candidates,
                      MyTracksData const& tracks)
  {
    for (const auto& candidate : candidates) {
      if (yCandMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
        continue;
      }
      if (candidate.pt() > ptCandMax) {
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
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFV0M()));
      int poolBinDs = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFV0M()));
      registry.fill(HIST("hTracksPoolBin"), poolBin);
      registry.fill(HIST("hDsPoolBin"), poolBinDs);
      for (const auto& [candidate, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (yCandMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandMax) {
          continue;
        }
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
          entryDsHadronGenInfo(isDsPrompt);
        } else if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          entryDsHadronPair(getDeltaPhi(pAssoc.phi(), candidate.phi()),
                            pAssoc.eta() - candidate.eta(),
                            candidate.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(hfHelper.invMassDsToPiKK(candidate), isDsSignal);
          entryDsHadronGenInfo(isDsPrompt);
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcRecME, "Process Mixed Event MCRec", false);

  // Event Mixing for the MCGen Mode (To do later)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDsHadronsSelCollision>(cfgc),
                      adaptAnalysisTask<HfCorrelatorDsHadrons>(cfgc)};
}
