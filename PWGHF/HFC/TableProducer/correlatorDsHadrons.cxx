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
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::aod::hf_correlation_ds_hadron;
using namespace o2::analysis::hf_cuts_ds_to_k_k_pi;
using namespace o2::constants::math;

namespace o2::aod
{
namespace ds_sel_collision
{
DECLARE_SOA_COLUMN(DsFound, dsFound, bool);
} // namespace hash
DECLARE_SOA_TABLE(DsSelCollision, "AOD", "DSCOLL", ds_sel_collision::DsFound);
using It = DsSelCollision::iterator;
} // namespace o2::aod

using namespace o2::aod::ds_sel_collision;

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

// binning type
std::vector<double> zBins{VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0};
std::vector<double> multBins{VARIABLE_WIDTH, 0., 200., 500.0, 5000.};
std::vector<double> multBinsMcGen{VARIABLE_WIDTH, 0., 20., 50.0, 500.}; // In MCGen multiplicity is defined by counting primaries
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
BinningType corrBinning{{zBins, multBins}, true};

struct HfDsSelectionCollision {
  Produces<aod::DsSelCollision> collisionSelDs; 

  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  Filter dsFlagFilter = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << DecayType::DsToKKPi)) != static_cast<uint8_t>(0); // filter in HfCand3Prong

  using candDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using candDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using candDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Partition<candDsData> selectedDsAllCand = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>> recoFlagDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;
     
  void processDsSelCollisionsData(aod::Collision const& collision, candDsData const& candidates)
  {
    bool isDsFound = false;
    if (selectedDsAllCand.size() > 0) {
      auto selectedDsAllCandGrouped = selectedDsAllCand->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      for (auto const& candidate : selectedDsAllCandGrouped) {
        if (yCandMax >= 0. && std::abs(yDplus(candidate)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
          continue;
        }
        isDsFound = true;
        break;
      }
    }
    collisionSelDs(isDsFound);
  }
  PROCESS_SWITCH(HfDsSelectionCollision, processDsSelCollisionsData, "Process Ds Collision Selection Data", true);

 void processDsSelCollisionsMcRec(aod::Collision const& collision, soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec> const& candidates)
  {
    bool isDsFound = false;
    if (recoFlagDsCandidates.size() > 0) {
      auto selectedDsCandidatesGroupedMc = recoFlagDsCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      for (auto const& candidate : selectedDsCandidatesGroupedMc) {
        // check decay channel flag for candidate
        if (!(candidate.hfflag() & 1 << DecayType::DplusToPiKPi)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(yDs(candidate)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
          continue;
        }
        isDsFound = true;
        break;
      }
    }
    collisionSelDs(isDsFound);
  }
  PROCESS_SWITCH(HfDsSelectionCollision, processDsSelCollisionsMcRec, "Process Ds Collision Selection MCRec", false);

  void processDsSelCollisionsMcGen(aod::McCollision const& mccollision, candDsMcGen const& particlesMc)
  {
    bool isDsFound = false;
    for (auto const& particle : particlesMc) {
      if (std::abs(particle.pdgCode()) != pdg::Code::kDS) {
        continue;
      }
      double yD = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
      if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
        continue;
      }
      isDsFound = true;
      break;
    }
    collisionSelDs(isDsFound);
  }
  PROCESS_SWITCH(HfDsSelectionCollision, processDsSelCollisionsMcGen, "Process Ds Collision Selection MCGen", false);
};

/// Ds-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorDsHadrons {
  SliceCache cache;
  Preslice<aod::HfCand3Prong> perCol = aod::hf_cand::collisionId;
  Produces<aod::DsHadronPair> entryDsHadronPair;
  Produces<aod::DsHadronRecoInfo> entryDsHadronRecoInfo;
  Produces<aod::DsHadronGenInfo> entryDsHadronGenInfo;

  Configurable<int> selectionFlagDs{"selectionFlagDs", 7, "Selection Flag for Ds"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCA_xy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCA_z of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 50., "max. track pT"};
  Configurable<float> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<float> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{vecEfficiencyDmeson}, "Efficiency values for Ds meson"};

  Filter collisionFilter = aod::ds_sel_collision::dsFound == true;
  using selCollisionsWithDs = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::DsSelCollision>>;
  using selCollisionsWithDsMc = soa::Join<aod::McCollisions, aod::Mults, aod::DsSelCollision>;

  Filter flagDsFilter = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << DecayType::DsToKKPi)) != static_cast<uint8_t>(0);
  using candDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using candDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using candDsMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Filter trackFilter = (aod::track::eta < std::abs(etaTrackMax)) && (aod::track::pt > ptTrackMin) 
                       && (aod::track::dcaXY < std::abs(dcaXYTrackMax)) && (aod::track::dcaZ < std::abs(dcaZTrackMax));
  using myTracksData = soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA>>;
  using myTracksMc = soa::Filtered<soa::Join<aod::BigTracksMC, aod::TracksDCA>>;
  
  Partition<candDsData> selectedDsAllCand = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;
  Partition<candDsData> selectedDsToKKPiCand = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs;
  Partition<candDsData> selectedDsToPiKKCand = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;
  Partition<candDsMcReco> selectedDsMcRecoCand = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "Ds,Hadron candidates", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng0", "Ds,Hadron candidates - prong 0", {HistType::kTH1F, {axisPtProng0}}},
     {"hPtProng1", "Ds,Hadron candidates - prong 1", {HistType::kTH1F, {axisPtProng1}}},
     {"hPtProng2", "Ds,Hadron candidates - prong 2", {HistType::kTH1F, {axisPtProng2}}},
     {"hSelectionStatusDsToKKPi", "Ds,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hSelectionStatusDsToPiKK", "Ds,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hEta", "Ds,Hadron candidates", {HistType::kTH1F, {axisEta}}},
     {"hPhi", "Ds,Hadron candidates", {HistType::kTH1F, {axisPhi}}},
     {"hY", "Ds,Hadron candidates", {HistType::kTH1F, {axisY}}},
     {"hPtCandMCRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtD}}},
     {"hPtCandMCRecSigPrompt", "Ds,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtD}}},
     {"hPtCandMCRecSigNonPrompt", "Ds,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng0MCRecSig", "Ds,Hadron candidates - MC Reco - prong 0", {HistType::kTH1F, {axisPtProng0}}},
     {"hPtProng1MCRecSig", "Ds,Hadron candidates - MC Reco - prong 1", {HistType::kTH1F, {axisPtProng1}}},
     {"hPtProng2MCRecSig", "Ds,Hadron candidates - MC Reco - prong 2", {HistType::kTH1F, {axisPtProng2}}},
     {"hPtCandMCRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtD}}},
     {"hPtProng0MCRecBkg", "Ds,Hadron candidates - MC Reco - prong 0", {HistType::kTH1F, {axisPtProng0}}},
     {"hPtProng1MCRecBkg", "Ds,Hadron candidates - MC Reco - prong 1", {HistType::kTH1F, {axisPtProng1}}},
     {"hPtProng2MCRecBkg", "Ds,Hadron candidates - MC Reco - prong 2", {HistType::kTH1F, {axisPtProng2}}},
     {"hSelectionStatusMCRec", "Ds,Hadron candidates - MC Reco;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hEtaMCRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}}},
     {"hPhiMCRecSig", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}}},
     {"hYMCRecSig", "Ds,Hadron candidates - MC Reco;;entries", {HistType::kTH1F, {axisY}}},
     {"hEtaMCRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}}},
     {"hPhiMCRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}}},
     {"hYMCRecBkg", "Ds,Hadron candidates - MC Reco", {HistType::kTH1F, {axisY}}},
     {"hMCEvtCount", "Event counter - MC Gen", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMCGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisPtD}}},
     {"hPtCandMCGenPrompt", "Ds,Hadron particles - MC Gen Prompt", {HistType::kTH1F, {axisPtD}}},
     {"hPtCandMCGenNonPrompt", "Ds,Hadron particles - MC Gen Non Prompt", {HistType::kTH1F, {axisPtD}}},
     {"hEtaMCGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisEta}}},
     {"hPhiMCGen", "Ds,Hadron particles - MC Gen", {HistType::kTH1F, {axisPhi}}},
     {"hYMCGen", "Ds,Hadron candidates - MC Gen", {HistType::kTH1F, {axisY}}},
     {"hcountDsHadronPerEvent", "Ds,Hadron particles - MC Gen;Number per event;entries", {HistType::kTH1F, {{21, -0.5, 20.5}}}},
     {"hMultiplicityPreSelection", "Multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hPtVsMultiplicityMCRecPrompt", "Multiplicity V0M - MC Rec Prompt", {HistType::kTH2F, {{axisPtD}, {1000, 0., 10000.}}}},
     {"hPtVsMultiplicityMCRecNonPrompt", "Multiplicity V0M - MC Rec Non Prompt", {HistType::kTH2F, {{axisPtD}, {1000, 0., 10000.}}}},
     {"hPtVsMultiplicityMCGenPrompt", "Multiplicity V0M - MC Gen Prompt", {HistType::kTH2F, {{axisPtD}, {1000, 0., 10000.}}}},
     {"hPtVsMultiplicityMCGenNonPrompt", "Multiplicity V0M - MC Gen Non Prompt", {HistType::kTH2F, {{axisPtD}, {1000, 0., 10000.}}}},
     {"hMultiplicity", "Multiplicity", {HistType::kTH1F, {{10000, 0., 10000., "multiplicity"}}}},
     {"hMultV0M", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hZvtx", "z vertex;z vertex;entries", {HistType::kTH1F, {{200, -20., 20.}}}},
     {"hDsPoolBin", "Ds selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}}},
     {"hTracksPoolBin", "Tracks selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}}}}};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMassDs2D", "Ds candidates", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDsData", "Ds candidates", {HistType::kTH1F, {axisMassD}});
    registry.add("hMassDsMCRec", "Ds candidates", {HistType::kTH1F, {axisMassD}});
    registry.add("hMassDsVsPtMCRec", "Ds signal candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDsMCRecSig", "Ds signal candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDsMCRecBkg", "Ds background candidates - MC Reco", {HistType::kTH2F, {{axisMassD}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hcountDstriggersMCGen", "Ds trigger particles - MC Gen", {HistType::kTH2F, {{1, -0.5, 0.5, "number of Ds triggers"}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
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
    registry.fill(HIST("hPhi"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf));
    registry.fill(HIST("hY"), yDs(candidate));
    return;
  }

  /// Fill histograms of quantities for the KKPi daugther-mass hypothesis for data
  /// \param candidate is candidate
  /// \param efficiencyWeight is the efficiency correction
  template <typename T1>
  void fillHistoKKPi(const T1& candidate, double efficiencyWeight)
  {
    registry.fill(HIST("hMassDs2D"), invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
    registry.fill(HIST("hMassDsData"), invMassDsToKKPi(candidate), efficiencyWeight);
    registry.fill(HIST("hSelectionStatusDsToKKPi"), candidate.isSelDsToKKPi());
    return;
  }

  /// Fill histograms of quantities for the PiKK daugther-mass hypothesis for data
  /// \param candidate is candidate
  /// \param efficiencyWeight is the efficiency correction
  template <typename T1>
  void fillHistoPiKK(const T1& candidate, double efficiencyWeight)
  {
    registry.fill(HIST("hMassDs2D"), invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
    registry.fill(HIST("hMassDsData"), invMassDsToPiKK(candidate), efficiencyWeight);
    registry.fill(HIST("hSelectionStatusDsToPiKK"), candidate.isSelDsToPiKK());
    return;
  }

  /// Fill histograms of quantities for the Ds signal for MC reco-level 
  /// \param candidate is candidate
  /// \param multiplicityV0M is the multiplicity
  template <typename T1>
  void fillHistoMcRecSig(const T1& candidate, float multiplicityV0M)
  {
    registry.fill(HIST("hPtCandMCRecSig"), candidate.pt());
    registry.fill(HIST("hPtProng0MCRecSig"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1MCRecSig"), candidate.ptProng1());
    registry.fill(HIST("hPtProng2MCRecSig"), candidate.ptProng2());
    registry.fill(HIST("hEtaMCRecSig"), candidate.eta());
    registry.fill(HIST("hPhiMCRecSig"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf));
    registry.fill(HIST("hYMCRecSig"), yDs(candidate));

    // prompt
    if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtCandMCRecSigPrompt"), candidate.pt());
      //std::cout <<"Multiplicity when Prompt: "<< multiplicityV0M << std::endl;
      registry.fill(HIST("hPtVsMultiplicityMCRecPrompt"), candidate.pt(), multiplicityV0M);
    }

    // non-prompt
    if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
      registry.fill(HIST("hPtCandMCRecSigNonPrompt"), candidate.pt());
      registry.fill(HIST("hPtVsMultiplicityMCRecNonPrompt"), candidate.pt(), multiplicityV0M);
    }

    return;
  }

  /// Fill histograms of quantities for the Ds backgroung for MC reco-level 
  /// \param candidate is candidate
  template <typename T1>
  void fillHistoMcRecBkg(const T1& candidate)
  {
    registry.fill(HIST("hPtCandMCRecBkg"), candidate.pt());
    registry.fill(HIST("hPtProng0MCRecBkg"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1MCRecBkg"), candidate.ptProng1());
    registry.fill(HIST("hPtProng2MCRecBkg"), candidate.ptProng2());
    registry.fill(HIST("hEtaMCRecBkg"), candidate.eta());
    registry.fill(HIST("hPhiMCRecBkg"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf));
    registry.fill(HIST("hYMCRecBkg"), yDs(candidate));
    return;
  }

  /// Fill histograms of quantities for the Ds signal for MC reco-level 
  /// \param particle is particle, Ds
  /// \param yD is the Ds rapidity
  template <typename T1>
  void fillHistoMcGen(const T1& particle, double yD)
  {
    registry.fill(HIST("hPtCandMCGen"), particle.pt());
    registry.fill(HIST("hEtaMCGen"), particle.eta());
    registry.fill(HIST("hPhiMCGen"), RecoDecay::constrainAngle(particle.phi(), -o2::constants::math::PIHalf));
    registry.fill(HIST("hYMCGen"), yD);
    registry.fill(HIST("hcountDstriggersMCGen"), 0, particle.pt()); // to count trigger Ds for normalisation

    // prompt
    if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
      registry.fill(HIST("hPtCandMCGenPrompt"), particle.pt());
    }

    // non-prompt
    if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
      registry.fill(HIST("hPtCandMCGenNonPrompt"), particle.pt());
    }

    return;
  }

  /// Ds-hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(selCollisionsWithDs::iterator const& collision, candDsData const& candidates, myTracksData const& tracks)
  {
    if (selectedDsAllCand.size() > 0) {
      registry.fill(HIST("hZvtx"), collision.posZ());
      registry.fill(HIST("hMultV0M"), collision.multFV0M());
      int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFV0M()));
      int nTracks = 0;
      if (collision.numContrib() > 1) {
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) { // se tolgo da errore perche track e' unused variable
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

      auto selectedDsToKKPiCandGrouped = selectedDsToKKPiCand->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      auto selectedDsToPiKKCandGrouped = selectedDsToPiKKCand->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      
      // Ds fill histograms and Ds-Hadron correlation for DsToKKPi
      for (const auto& candidate : selectedDsToKKPiCandGrouped) {
        if (yCandMax >= 0. && std::abs(yDs(candidate)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
          continue;
        }
        if (candidate.pt() > ptTrackMax) {
          continue;
        }
        double efficiencyWeight = 1.;
        if (applyEfficiency) {
          efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate.pt()));
        }
        fillHisto(candidate);
        fillHistoKKPi(candidate, efficiencyWeight);

        // Ds-Hadron correlation dedicated section
        for (const auto& track : tracks) {
          // Removing Ds daughters by checking track indices
          if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
            continue;
          }
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin);
          entryDsHadronRecoInfo(invMassDsToKKPi(candidate), false, false);
        }
      }

      // Ds fill histograms and Ds-Hadron correlation for DsToPiKK
      for (const auto& candidate : selectedDsToPiKKCandGrouped) {
        if (yCandMax >= 0. && std::abs(yDs(candidate)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
          continue;
        }
        if (candidate.pt() > ptTrackMax) {
          continue;
        }
        double efficiencyWeight = 1.;
        if (applyEfficiency) {
          efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate.pt()));
        }
        fillHisto(candidate);
        fillHistoPiKK(candidate, efficiencyWeight);

        // Ds-Hadron correlation dedicated section
        for (const auto& track : tracks) {
          // Removing Ds daughters by checking track indices
          if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
            continue;
          }
          entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin);
          entryDsHadronRecoInfo(invMassDsToPiKK(candidate), false, false);
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processData, "Process data", true);

  /// Ds-Hadron correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(selCollisionsWithDs::iterator const& collision, candDsMcReco const& candidates, candDsMcGen const& particlesMc, myTracksMc const& tracks)
  {
    if (selectedDsMcRecoCand.size() > 0) {
      registry.fill(HIST("hZvtx"), collision.posZ());
      registry.fill(HIST("hMultV0M"), collision.multFV0M());
      int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFV0M()));
      int nTracks = 0;
      if (collision.numContrib() > 1) {
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
            continue;
          }
          nTracks++;
        }
      }
      //registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
      if (nTracks < multMin || nTracks > multMax) {
        return;
      }
      registry.fill(HIST("hMultiplicity"), nTracks);

      auto selectedDsMcRecoCandGrouped = selectedDsMcRecoCand->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
       
      // MC reco level
      bool flagDsPrompt = false;
      float multiplicityV0M = collision.multFV0M();
      for (const auto& candidate : candidates) {
        if (yCandMax >= 0. && std::abs(yDs(candidate)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
          continue;
        }
        if (candidate.pt() >= ptTrackMax) {
          continue;
        }
        double efficiencyWeight = 1.;
        if (applyEfficiency) {
          efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate.pt()));
        }
        auto prong0McPart = candidate.prong0_as<myTracksMc>().mcParticle_as<candDsMcGen>();
        auto indexMother = RecoDecay::getMother(particlesMc, prong0McPart, pdg::Code::kDS, true);
        auto particleMother = particlesMc.iteratorAt(indexMother);
        if (std::abs(candidate.flagMcMatchRec()) == 1 << DecayType::DsToKKPi) {
          fillHistoMcRecSig(candidate, multiplicityV0M);
          // KKPi
          if (std::abs(prong0McPart.pdgCode()) == kKPlus) {
            registry.fill(HIST("hMassDsMCRec"), invMassDsToKKPi(candidate), efficiencyWeight);
            registry.fill(HIST("hMassDsMCRecSig"), invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hMassDsVsPtMCRec"), invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hSelectionStatusMCRec"), candidate.isSelDsToKKPi());
          }
          // PiKK
          if (std::abs(prong0McPart.pdgCode()) == kPiPlus) {
            registry.fill(HIST("hMassDsMCRec"), invMassDsToPiKK(candidate), efficiencyWeight);
            registry.fill(HIST("hMassDsMCRecSig"), invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hMassDsVsPtMCRec"), invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hSelectionStatusMCRec"), candidate.isSelDsToPiKK());
          }
          // prompt
          if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          flagDsPrompt = true;
          }
          // non-prompt
          if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
          flagDsPrompt = false;
          }
        } else {
          fillHistoMcRecBkg(candidate);
          // KKPi
          if (std::abs(prong0McPart.pdgCode()) == kKPlus) {
            registry.fill(HIST("hMassDsMCRec"), invMassDsToKKPi(candidate), efficiencyWeight);
            registry.fill(HIST("hMassDsMCRecBkg"), invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hMassDsVsPtMCRec"), invMassDsToKKPi(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hSelectionStatusMCRec"), candidate.isSelDsToKKPi());
          }
          // PiKK
          if (std::abs(prong0McPart.pdgCode()) == kPiPlus) {
            registry.fill(HIST("hMassDsMCRec"), invMassDsToPiKK(candidate), efficiencyWeight);
            registry.fill(HIST("hMassDsMCRecBkg"), invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hMassDsVsPtMCRec"), invMassDsToPiKK(candidate), candidate.pt(), efficiencyWeight);
            registry.fill(HIST("hSelectionStatusMCRec"), candidate.isSelDsToPiKK());
          }
        }

        // Ds-Hadron correlation dedicated section
        // if the candidate is selected as Ds, search for Hadron and evaluate correlations
        bool flagDsSignal = std::abs(candidate.flagMcMatchRec()) == 1 << DecayType::DsToKKPi;
        for (const auto& track : tracks) {
          // Removing Ds daughters by checking track indices
          if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
            continue;
          }
          // KKPi
          if (std::abs(prong0McPart.pdgCode()) == kKPlus) {
            entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                              track.eta() - candidate.eta(),
                              candidate.pt(),
                              track.pt(),
                              poolBin);
            entryDsHadronRecoInfo(invMassDsToKKPi(candidate), flagDsSignal, flagDsPrompt);
          }
          // PiKK
          if (std::abs(prong0McPart.pdgCode()) == kPiPlus) {
            entryDsHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                              track.eta() - candidate.eta(),
                              candidate.pt(),
                              track.pt(),
                              poolBin);
            entryDsHadronRecoInfo(invMassDsToPiKK(candidate), flagDsSignal, flagDsPrompt);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcRec, "Process MC Reco mode", false);

  /// Ds-Hadron correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(selCollisionsWithDsMc::iterator const& mccollision, candDsMcGen const& particlesMc)
  {
    int counterDsHadron = 0;
    registry.fill(HIST("hMCEvtCount"), 0);

    auto getTracksSize = [&particlesMc](selCollisionsWithDsMc::iterator const& mccollision) {
      int nTracks = 0;
      for (auto& track : particlesMc) {
        if (track.isPhysicalPrimary() && std::abs(track.eta()) < 1.0) {
          nTracks++;
        }
      }
      return nTracks;
    };
    using BinningTypeMcGen = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::mccollision::PosZ, decltype(getTracksSize)>;
    BinningTypeMcGen corrBinningMcGen{{getTracksSize}, {zBins, multBinsMcGen}, true};
    int poolBin = corrBinningMcGen.getBin(std::make_tuple(mccollision.posZ(), getTracksSize(mccollision)));

    // MC gen level
    for (auto const& particle : particlesMc) {
      // check if the particle is Ds
      if (std::abs(particle.pdgCode()) != pdg::Code::kDS) {
        continue;
      }
      if (std::abs(particle.flagMcMatchGen()) == 1 << DecayType::DsToKKPi){ 
        double yD = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
          continue;
        }
        fillHistoMcGen(particle, yD);
        counterDsHadron++;
        bool flagDsPrompt = false;
        // prompt
        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          flagDsPrompt = true;
        }
        // non-prompt
        if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          flagDsPrompt = false;
        }
        // Ds Hadron correlation dedicated section
        for (auto const& particleAssoc : particlesMc) {
          if (std::abs(particleAssoc.eta()) > etaTrackMax) {
            continue;
          }
          if (particleAssoc.pt() < ptTrackMin) {
            continue;
          }

          if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particleAssoc.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
            continue;
          }

          entryDsHadronPair(getDeltaPhi(particleAssoc.phi(), particle.phi()),
                            particleAssoc.eta() - particle.eta(),
                            particle.pt(),
                            particleAssoc.pt(),
                            poolBin);
          entryDsHadronGenInfo(flagDsPrompt);
        }
      }
    }
    registry.fill(HIST("hcountDsHadronPerEvent"), counterDsHadron);
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcGen, "Process MC Gen mode", false);

  // Event Mixing
   void processDataME(selCollisionsWithDs& collisions, candDsData& candidates, myTracksData& tracks)
    {
      auto tracksTuple = std::make_tuple(candidates,tracks);
      Pair<selCollisionsWithDs, candDsData, myTracksData, BinningType> pairData{corrBinning, 5, -1, collisions, tracksTuple, &cache};
   
    for (auto& [c1, tracks1, c2, tracks2] : pairData) {
       LOGF(info, "Mixed event collisions: Index = (%d, %d), tracks Size: (%d, %d), Z Vertex: (%f, %f), Pool Bin: (%d, %d)", c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size(), c1.posZ(), c2.posZ(), corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFV0M())),corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFV0M())));
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFV0M()));
      for (auto& [candidate, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) { 
          
        if (yCandMax >= 0. && std::abs(yDs(candidate)) > yCandMax) {
          continue;
        }
        std::cout << "candidate.isSelDsToKKPi() = " << candidate.isSelDsToKKPi() << std::endl;
        std::cout << "candidate.isSelDsToPiKK() = " << candidate.isSelDsToPiKK() << std::endl;
        // KKPi
        if (candidate.isSelDsToKKPi() == selectionFlagDs){
          LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", candidate.index(), pAssoc.index(), c1.index(), c2.index(), candidate.collision().index(), pAssoc.collision().index()); 
          entryDsHadronPair(getDeltaPhi(candidate.phi(), pAssoc.phi()), 
                            candidate.eta() - pAssoc.eta(), 
                            candidate.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(invMassDsToKKPi(candidate), false, false);
        }
        // PiKK
        else if (candidate.isSelDsToPiKK() == selectionFlagDs){
          LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", candidate.index(), pAssoc.index(), c1.index(), c2.index(), candidate.collision().index(), pAssoc.collision().index()); 
          entryDsHadronPair(getDeltaPhi(candidate.phi(), pAssoc.phi()), 
                            candidate.eta() - pAssoc.eta(), 
                            candidate.pt(),
                            pAssoc.pt(),
                            poolBin);
          entryDsHadronRecoInfo(invMassDsToPiKK(candidate), false, false);
        }
        else {
          LOGF(info, "**** PROBLEM IN THE DS SELECTION OR IN THE FLAG USED ****");
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processDataME, "Process Mixed Event Data", false);

  void processMcRecME(selCollisionsWithDs& collisions, candDsMcReco& candidates, candDsMcGen const& particlesMc, myTracksMc& tracks)
  {
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<selCollisionsWithDs, candDsMcReco, myTracksMc, BinningType> pairMcRec{corrBinning, 5, -1, collisions, tracksTuple, &cache};
    
    //bool flagDsSignal = false;
    bool flagDsPrompt = false;
    for (auto& [c1, tracks1, c2, tracks2] : pairMcRec) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFV0M()));
      for (auto& [candidate, pAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (yCandMax >= 0. && std::abs(yDplus(candidate)) > yCandMax) {
          continue;
        }
        auto prong0McPart = candidate.prong0_as<myTracksMc>().mcParticle_as<candDsMcGen>();
        bool flagDsSignal = std::abs(candidate.flagMcMatchRec()) == 1 << DecayType::DsToKKPi;
        // prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          flagDsPrompt = true;
        }
        // non-prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
          flagDsPrompt = false;
        }
        // KKPi
        if (std::abs(prong0McPart.pdgCode()) == kKPlus) {
          entryDsHadronPair(getDeltaPhi(candidate.phi(), pAssoc.phi()),
                          candidate.eta() - pAssoc.eta(),
                          candidate.pt(),
                          pAssoc.pt(),
                          poolBin);
          entryDsHadronRecoInfo(invMassDsToKKPi(candidate), flagDsSignal, flagDsPrompt);
        }
        // PiKK
        if (std::abs(prong0McPart.pdgCode()) == kPiPlus) {
          entryDsHadronPair(getDeltaPhi(candidate.phi(), pAssoc.phi()),
                          candidate.eta() - pAssoc.eta(),
                          candidate.pt(),
                          pAssoc.pt(),
                          poolBin);
          entryDsHadronRecoInfo(invMassDsToPiKK(candidate), flagDsSignal, flagDsPrompt);
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDsHadrons, processMcRecME, "Process Mixed Event MCRec", false);

  // Event Mixing for the MCGen Mode (To do later)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDsSelectionCollision>(cfgc),
                      adaptAnalysisTask<HfCorrelatorDsHadrons>(cfgc)};
}
