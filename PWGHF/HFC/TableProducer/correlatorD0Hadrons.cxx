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

/// \file correlatorD0Hadrons.cxx
/// \brief D0-Hadron correlator task - data-like, MC-reco and MC-kine analyses.
///
/// \author Samrangy Sadhu <samrangy.sadhu@cern.ch>, INFN Bari
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>, IIT Indore

#include <vector>

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
#include "PWGHF/HFC/Utils/utilsCorrelations.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::hf_correlations;

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiHadron, double phiD)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PIHalf);
}

const int nPtBinsMassAndEfficiency = o2::analysis::hf_cuts_d0_to_pi_k::nBinsPt;
const double efficiencyDmesonDefault[nPtBinsMassAndEfficiency] = {};
auto vecEfficiencyDmeson = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + nPtBinsMassAndEfficiency};

// histogram binning definition
const int massAxisNBins = 200;
const double massAxisMin = 1.3848;
const double massAxisMax = 2.3848;
const int phiAxisNBins = 32;
const double phiAxisMin = 0.;
const double phiAxisMax = o2::constants::math::TwoPI;
const int yAxisNBins = 100;
const double yAxisMin = -5.;
const double yAxisMax = 5.;
const int ptDAxisNBins = 180;
const double ptDAxisMin = 0.;
const double ptDAxisMax = 36.;

// Types
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::DmesonSelection>>;
using SelectedTracks = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>;
using SelectedCandidatesData = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
using SelectedTracksMcRec = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels>>;
using SelectedCandidatesMcRec = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>;
using SelectedCollisionsMcGen = soa::Filtered<soa::Join<aod::McCollisions, aod::DmesonSelection, aod::MultsExtraMC>>;
using SelectedParticlesMcGen = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;

// Code to select collisions with at least one D0
struct HfCorrelatorD0HadronsSelection {
  SliceCache cache;

  Produces<aod::DmesonSelection> d0Sel;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<float> yCandMax{"yCandMax", 4.0, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", -1., "min. cand. pT"};

  HfHelper hfHelper;

  Preslice<aod::HfCand2Prong> perCol = aod::hf_cand::collisionId;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> selectedD0candidatesMc = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  void processD0SelectionData(aod::Collision const& collision,
                              soa::Join<aod::HfCand2Prong, aod::HfSelD0> const&)
  {
    bool isD0Found = 0;
    if (selectedD0Candidates.size() > 0) {
      auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

      for (const auto& candidate : selectedD0CandidatesGrouped) {
        // check decay channel flag for candidate
        if (!TESTBIT(candidate.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          continue;
        }
        if (std::abs(hfHelper.yD0(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          continue;
        }
        isD0Found = 1;
        break;
      }
    }
    d0Sel(isD0Found);
  }
  PROCESS_SWITCH(HfCorrelatorD0HadronsSelection, processD0SelectionData, "Process D0 Selection Data", false);

  void processD0SelectionMcRec(aod::Collision const& collision,
                               soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const&)
  {
    bool isD0Found = 0;
    if (selectedD0candidatesMc.size() > 0) {
      auto selectedD0CandidatesGroupedMc = selectedD0candidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      for (const auto& candidate : selectedD0CandidatesGroupedMc) {
        // check decay channel flag for candidate
        if (!TESTBIT(candidate.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          continue;
        }
        if (std::abs(hfHelper.yD0(candidate)) > yCandMax || candidate.pt() < ptCandMin) {
          continue;
        }
        isD0Found = 1;
        break;
      }
    }
    d0Sel(isD0Found);
  }
  PROCESS_SWITCH(HfCorrelatorD0HadronsSelection, processD0SelectionMcRec, "Process D0 Selection MCRec", true);

  void processD0SelectionMcGen(aod::McCollision const&,
                               aod::McParticles const& mcParticles)
  {
    bool isD0Found = 0;
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != Pdg::kD0) {
        continue;
      }
      double yD = RecoDecay::y(particle.pVector(), MassD0);
      if (std::abs(yD) > yCandMax || particle.pt() < ptCandMin) {
        continue;
      }
      isD0Found = 1;
      break;
    }
    d0Sel(isD0Found);
  }
  PROCESS_SWITCH(HfCorrelatorD0HadronsSelection, processD0SelectionMcGen, "Process D0 Selection MCGen", false);
};

/// D0-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorD0Hadrons {
  SliceCache cache;

  Produces<aod::DHadronPair> entryD0HadronPair;
  Produces<aod::DHadronRecoInfo> entryD0HadronRecoInfo;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCAxy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCAz of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 99., "max. track pT"};
  Configurable<std::vector<double>> bins{"ptBinsForMassAndEfficiency", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyDmeson{"efficiencyDmeson", std::vector<double>{vecEfficiencyDmeson}, "Efficiency values for D0 meson"};
  Configurable<int> applyEfficiency{"efficiencyFlagD", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<float> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<float> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<float> ptSoftPionMax{"ptSoftPionMax", 3.f * 800.f * std::pow(10.f, -6.f), "max. pT cut for soft pion identification"};
  Configurable<bool> correlateD0WithLeadingParticle{"correlateD0WithLeadingParticle", false, "Switch for correlation of D0 mesons with leading particle only"};
  Configurable<bool> storeAutoCorrelationFlag{"storeAutoCorrelationFlag", false, "Store flag that indicates if the track is paired to its D-meson mother instead of skipping it"};
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "z vertex position pools"};
  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 10000.0f}, "event multiplicity pools (FT0M)"};
  ConfigurableAxis multPoolBinsMcGen{"multPoolBinsMcGen", {VARIABLE_WIDTH, 0.0f, 20.0f, 50.0f, 500.0f}, "Mixing bins - MC multiplicity"}; // In MCGen multiplicity is defined by counting tracks

  HfHelper hfHelper;
  BinningType corrBinning{{zPoolBins, multPoolBins}, true};

  int leadingIndex = 0;
  double massD0{0.};
  double massPi{0.};
  double massK{0.};
  double softPiMass = 0.14543; // pion mass + Q-value of the D*->D0pi decay

  Preslice<aod::HfCand2Prong> perCol = aod::hf_cand::collisionId;

  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> selectedD0candidatesMc = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter trackFilter = requireGlobalTrackWoDCAInFilter() && (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);
  Filter d0Filter = (aod::hf_sel_candidate_d0::isSelD0 >= 1) || (aod::hf_sel_candidate_d0::isSelD0bar >= 1);
  Filter collisionFilterGen = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter particlesFilter = nabs(aod::mcparticle::pdgCode) == static_cast<int>(Pdg::kD0) || ((aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary);

  HistogramRegistry registry{
    "registry",
    // NOTE: use hMassD0 for trigger normalisation (S*0.955), and hMass2DCorrelationPairs (in final task) for 2D-sideband-subtraction purposes
    {{"hPtCand", "D0,D0bar candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0", "D0,D0bar candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1", "D0,D0bar candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatus", "D0,D0bar candidates;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEta", "D0,D0bar candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hPhi", "D0,D0bar candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisNBins, phiAxisMin, phiAxisMax}}}},
     {"hY", "D0,D0bar candidates;candidate #it{y};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hPtCandRec", "D0,D0bar candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0Rec", "D0,D0bar candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1Rec", "D0,D0bar candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatusRec", "D0,D0bar candidates - MC reco;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hSignalStatusMERec", "Signal Status - MC reco ME;candidate sidnalStatus;entries", {HistType::kTH1F, {{200, 0, 200}}}},
     {"hEtaRec", "D0,D0bar candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hPhiRec", "D0,D0bar candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisNBins, phiAxisMin, phiAxisMax}}}},
     {"hYRec", "D0,D0bar candidates - MC reco;candidate #it{y};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hEvtCountGen", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandGen", "D0,D0bar particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisNBins, ptDAxisMin, ptDAxisMax}}}},
     {"hEtaGen", "D0,D0bar particles - MC gen;particle #it{#eta};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hPhiGen", "D0,D0bar particles - MC gen;particle #it{#varphi};entries", {HistType::kTH1F, {{phiAxisNBins, phiAxisMin, phiAxisMax}}}},
     {"hYGen", "D0,D0bar candidates - MC gen;candidate #it{y};entries", {HistType::kTH1F, {{yAxisNBins, yAxisMin, yAxisMax}}}},
     {"hTrackCounter", "soft pion counter -  Data", {HistType::kTH1F, {{5, 0., 5.}}}},
     {"hTrackCounterRec", "soft pion counter - MC rec", {HistType::kTH1F, {{5, 0., 5.}}}},
     {"hTrackCounterGen", "soft pion counter - MC gen", {HistType::kTH1F, {{5, 0., 5.}}}},
     {"hMultV0M", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hZvtx", "z vertex;z vertex;entries", {HistType::kTH1F, {{200, -20., 20.}}}},
     {"hMultFT0M", "Multiplicity FT0M", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hD0Bin", "D0 selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}}},
     {"hTracksBin", "Tracks selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}}}}};

  void init(InitContext&)
  {
    massD0 = MassD0;
    massPi = MassPiPlus;
    massK = MassKPlus;

    auto vbins = (std::vector<double>)bins;
    registry.add("hMass", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMass1D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisNBins, massAxisMin, massAxisMax}}});
    registry.add("hMassD01D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisNBins, massAxisMin, massAxisMax}}});
    registry.add("hMassD0bar1D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisNBins, massAxisMin, massAxisMax}}});
    // mass histogram for D0 signal candidates
    registry.add("hMassD0RecSig", "D0 signal candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for D0 Reflection candidates
    registry.add("hMassD0RecRef", "D0 reflection candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for D0 background candidates
    registry.add("hMassD0RecBg", "D0 background candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for D0bar signal candidates
    registry.add("hMassD0barRecSig", "D0bar signal candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for D0bar Reflection candidates
    registry.add("hMassD0barRecRef", "D0bar reflection candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for D0bar background candidates
    registry.add("hMassD0barRecBg", "D0bar background candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisNBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCountD0TriggersGen", "D0 trigger particles - MC gen;;N of trigger D0", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  // =======  Process starts for Data, Same event ============

  /// D0-h correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(SelectedCollisions::iterator const& collision,
                   SelectedTracks const& tracks,
                   soa::Join<aod::HfCand2Prong, aod::HfSelD0> const&)
  {
    // protection against empty tables to be sliced
    if (selectedD0Candidates.size() == 0) {
      return;
    }
    // find leading particle
    if (correlateD0WithLeadingParticle) {
      leadingIndex = findLeadingParticle(tracks, dcaXYTrackMax.value, dcaZTrackMax.value, etaTrackMax.value);
    }

    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
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
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), collision.multFT0M());
    if (nTracks < multMin || nTracks > multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

    for (const auto& candidate1 : selectedD0CandidatesGrouped) {
      if (std::abs(hfHelper.yD0(candidate1)) >= yCandMax || candidate1.pt() <= ptCandMin || candidate1.pt() >= ptTrackMax) {
        continue;
      }
      // check decay channel flag for candidate1
      if (!TESTBIT(candidate1.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }

      // ========================== Define parameters for soft pion removal ================================
      auto ePiK = RecoDecay::e(candidate1.pVectorProng0(), massPi) + RecoDecay::e(candidate1.pVectorProng1(), massK);
      auto eKPi = RecoDecay::e(candidate1.pVectorProng0(), massK) + RecoDecay::e(candidate1.pVectorProng1(), massPi);

      // ========================== trigger efficiency ================================
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(bins, candidate1.pt()));
      }
      // ========================== Fill mass histo  ================================
      if (candidate1.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMass1D"), hfHelper.invMassD0ToPiK(candidate1), efficiencyWeight);
        registry.fill(HIST("hMassD01D"), hfHelper.invMassD0ToPiK(candidate1), efficiencyWeight);
      }
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMass1D"), hfHelper.invMassD0barToKPi(candidate1), efficiencyWeight);
        registry.fill(HIST("hMassD0bar1D"), hfHelper.invMassD0barToKPi(candidate1), efficiencyWeight);
      }
      // ========================== Fill general histos ================================
      registry.fill(HIST("hPtCand"), candidate1.pt());
      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), hfHelper.yD0(candidate1));
      registry.fill(HIST("hSelectionStatus"), candidate1.isSelD0bar() + (candidate1.isSelD0() * 2));
      registry.fill(HIST("hD0Bin"), poolBin);

      // ============ D-h correlation dedicated section ==================================

      // ========================== track loop starts here ================================
      for (const auto& track : tracks) {
        registry.fill(HIST("hTrackCounter"), 1); // fill total no. of tracks
        // Remove D0 daughters by checking track indices
        bool correlationStatus = false;
        if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex())) {
          if (!storeAutoCorrelationFlag) {
            continue;
          }
          correlationStatus = true;
        }

        registry.fill(HIST("hTrackCounter"), 2); // fill no. of tracks before soft pion removal

        // ========== soft pion removal ===================================================
        double invMassDstar1 = 0., invMassDstar2 = 0.;
        bool isSoftPiD0 = false, isSoftPiD0bar = false;
        auto pSum2 = RecoDecay::p2(candidate1.pVector(), track.pVector());
        auto ePion = track.energy(massPi);
        invMassDstar1 = std::sqrt((ePiK + ePion) * (ePiK + ePion) - pSum2);
        invMassDstar2 = std::sqrt((eKPi + ePion) * (eKPi + ePion) - pSum2);

        if (candidate1.isSelD0() >= selectionFlagD0) {
          if ((std::abs(invMassDstar1 - hfHelper.invMassD0ToPiK(candidate1)) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0 = true;
            continue;
          }
        }

        if (candidate1.isSelD0bar() >= selectionFlagD0bar) {
          if ((std::abs(invMassDstar2 - hfHelper.invMassD0barToKPi(candidate1)) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0bar = true;
            continue;
          }
        }
        registry.fill(HIST("hTrackCounter"), 3); // fill no. of tracks after soft pion removal

        int signalStatus = 0;
        if ((candidate1.isSelD0() >= selectionFlagD0) && !isSoftPiD0) {
          signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0Only;
        }
        if ((candidate1.isSelD0bar() >= selectionFlagD0bar) && !isSoftPiD0bar) {
          signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0barOnly;
        }

        if (correlateD0WithLeadingParticle) {
          if (track.globalIndex() != leadingIndex) {
            continue;
          }
          registry.fill(HIST("hTrackCounter"), 4); // fill no. of tracks  have leading particle
        }
        entryD0HadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                          track.eta() - candidate1.eta(),
                          candidate1.pt(),
                          track.pt(),
                          poolBin,
                          correlationStatus);
        entryD0HadronRecoInfo(hfHelper.invMassD0ToPiK(candidate1), hfHelper.invMassD0barToKPi(candidate1), signalStatus);

      } // end inner loop (tracks)

    } // end outer loop
  }
  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processData, "Process data", false);

  // ================  Process starts for MCRec, same event ========================

  void processMcRec(SelectedCollisions::iterator const& collision,
                    SelectedTracksMcRec const& tracks,
                    soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const&)
  {
    // protection against empty tables to be sliced
    if (selectedD0candidatesMc.size() == 0) {
      return;
    }
    // find leading particle
    if (correlateD0WithLeadingParticle) {
      leadingIndex = findLeadingParticle(tracks, dcaXYTrackMax.value, dcaZTrackMax.value, etaTrackMax.value);
    }
    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
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
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < multMin || nTracks > multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    auto selectedD0CandidatesGroupedMc = selectedD0candidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    // MC reco level
    bool flagD0 = false;
    bool flagD0bar = false;

    for (const auto& candidate1 : selectedD0CandidatesGroupedMc) {
      // check decay channel flag for candidate1
      if (!TESTBIT(candidate1.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (std::abs(hfHelper.yD0(candidate1)) >= yCandMax || candidate1.pt() <= ptCandMin || candidate1.pt() >= ptTrackMax) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(bins, candidate1.pt()));
      }

      if (std::abs(candidate1.flagMcMatchRec()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) {
        // fill per-candidate distributions from D0/D0bar true candidates
        registry.fill(HIST("hPtCandRec"), candidate1.pt());
        registry.fill(HIST("hPtProng0Rec"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1Rec"), candidate1.ptProng1());
        registry.fill(HIST("hEtaRec"), candidate1.eta());
        registry.fill(HIST("hPhiRec"), candidate1.phi());
        registry.fill(HIST("hYRec"), hfHelper.yD0(candidate1));
        registry.fill(HIST("hSelectionStatusRec"), candidate1.isSelD0bar() + (candidate1.isSelD0() * 2));
      }
      // fill invariant mass plots from D0/D0bar signal and background candidates
      if (candidate1.isSelD0() >= selectionFlagD0) {                                       // only reco as D0
        if (candidate1.flagMcMatchRec() == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) { // also matched as D0
          registry.fill(HIST("hMassD0RecSig"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        } else if (candidate1.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassD0RecRef"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0RecBg"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt(), efficiencyWeight);
        }
      }
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) {                                    // only reco as D0bar
        if (candidate1.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) { // also matched as D0bar
          registry.fill(HIST("hMassD0barRecSig"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else if (candidate1.flagMcMatchRec() == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) {
          registry.fill(HIST("hMassD0barRecRef"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0barRecBg"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt(), efficiencyWeight);
        }
      }

      // ===================== Define parameters for soft pion removal ========================
      auto ePiK = RecoDecay::e(candidate1.pVectorProng0(), massPi) + RecoDecay::e(candidate1.pVectorProng1(), massK);
      auto eKPi = RecoDecay::e(candidate1.pVectorProng0(), massK) + RecoDecay::e(candidate1.pVectorProng1(), massPi);

      // ============== D-h correlation dedicated section ====================================

      flagD0 = candidate1.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK);     // flagD0Signal 'true' if candidate1 matched to D0 (particle)
      flagD0bar = candidate1.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK); // flagD0Reflection 'true' if candidate1, selected as D0 (particle), is matched to D0bar (antiparticle)

      // ========== track loop starts here ========================

      for (const auto& track : tracks) {
        registry.fill(HIST("hTrackCounterRec"), 1); // fill total no. of tracks

        // Removing D0 daughters by checking track indices
        bool correlationStatus = false;
        if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex())) {
          if (!storeAutoCorrelationFlag) {
            continue;
          }
          correlationStatus = true;
        }
        registry.fill(HIST("hTrackCounterRec"), 2); // fill no. of tracks before soft pion removal

        // ===== soft pion removal ===================================================
        double invMassDstar1 = 0, invMassDstar2 = 0;
        bool isSoftPiD0 = false, isSoftPiD0bar = false;
        auto pSum2 = RecoDecay::p2(candidate1.pVector(), track.pVector());
        auto ePion = track.energy(massPi);
        invMassDstar1 = std::sqrt((ePiK + ePion) * (ePiK + ePion) - pSum2);
        invMassDstar2 = std::sqrt((eKPi + ePion) * (eKPi + ePion) - pSum2);

        if (candidate1.isSelD0() >= selectionFlagD0) {
          if ((std::abs(invMassDstar1 - hfHelper.invMassD0ToPiK(candidate1)) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0 = true;
            continue;
          }
        }

        if (candidate1.isSelD0bar() >= selectionFlagD0bar) {
          if ((std::abs(invMassDstar2 - hfHelper.invMassD0barToKPi(candidate1)) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0bar = true;
            continue;
          }
        }

        registry.fill(HIST("hTrackCounterRec"), 3); // fill no. of tracks after soft pion removal

        if (correlateD0WithLeadingParticle) {
          if (track.globalIndex() != leadingIndex) {
            continue;
          }
          registry.fill(HIST("hTrackCounterRec"), 4); // fill no. of tracks  have leading particle
        }

        int signalStatus = 0;
        if (flagD0 && (candidate1.isSelD0() >= selectionFlagD0) && !isSoftPiD0) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Sig);
        } // signal case D0
        if (flagD0bar && (candidate1.isSelD0() >= selectionFlagD0) && !isSoftPiD0) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Ref);
        } // reflection case D0
        if (!flagD0 && !flagD0bar && (candidate1.isSelD0() >= selectionFlagD0) && !isSoftPiD0) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Bg);
        } // background case D0

        if (flagD0bar && (candidate1.isSelD0bar() >= selectionFlagD0bar) && !isSoftPiD0bar) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barSig);
        } // signal case D0bar
        if (flagD0 && (candidate1.isSelD0bar() >= selectionFlagD0bar) && !isSoftPiD0bar) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barRef);
        } // reflection case D0bar
        if (!flagD0 && !flagD0bar && (candidate1.isSelD0bar() >= selectionFlagD0bar) && !isSoftPiD0bar) {
          SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barBg);
        } // background case D0bar

        entryD0HadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                          track.eta() - candidate1.eta(),
                          candidate1.pt(),
                          track.pt(),
                          poolBin,
                          correlationStatus);
        entryD0HadronRecoInfo(hfHelper.invMassD0ToPiK(candidate1), hfHelper.invMassD0barToKPi(candidate1), signalStatus);
      } // end inner loop (Tracks)
    }   // end of outer loop (D0)
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultV0M"), collision.multFT0M());
  }

  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processMcRec, "Process MC Reco mode", true);

  // =================  Process starts for MCGen, same event ===================

  void processMcGen(aod::McCollision const& mcCollision,
                    soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles)
  {
    registry.fill(HIST("hEvtCountGen"), 0);
    // MC gen level
    // find leading particle
    if (correlateD0WithLeadingParticle) {
      leadingIndex = findLeadingParticleMcGen(mcParticles, etaTrackMax.value, ptTrackMin.value);
    }
    for (const auto& particle1 : mcParticles) {
      // check if the particle is D0 or D0bar (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != Pdg::kD0) {
        continue;
      }
      double yD = RecoDecay::y(particle1.pVector(), MassD0);
      if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle1.pt() < ptCandMin) {
        continue;
      }
      registry.fill(HIST("hPtCandGen"), particle1.pt());
      registry.fill(HIST("hEtaGen"), particle1.eta());
      registry.fill(HIST("hPhiGen"), particle1.phi());
      registry.fill(HIST("hYGen"), yD);

      // =============== D-h correlation dedicated section =====================

      if (std::abs(particle1.pdgCode()) != Pdg::kD0) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
        continue;
      }
      registry.fill(HIST("hCountD0TriggersGen"), 0, particle1.pt()); // to count trigger D0 (for normalisation)

      for (const auto& particle2 : mcParticles) {
        registry.fill(HIST("hTrackCounterGen"), 1); // total no. of tracks
        if (std::abs(particle2.eta()) > etaTrackMax) {
          continue;
        }
        if (particle2.pt() < ptTrackMin) {
          continue;
        }
        if ((std::abs(particle2.pdgCode()) != kElectron) && (std::abs(particle2.pdgCode()) != kMuonMinus) && (std::abs(particle2.pdgCode()) != kPiPlus) && (std::abs(particle2.pdgCode()) != kKPlus) && (std::abs(particle2.pdgCode()) != kProton)) {
          continue;
        }

        // ==============================soft pion removal================================
        registry.fill(HIST("hTrackCounterGen"), 2); // fill before soft pi removal
        // method used: indexMother = -1 by default if the mother doesn't match with given PID of the mother. We find mother of pion if it is D* and mother of D0 if it is D*. If they are both positive and they both match each other, then it is detected as a soft pion

        auto indexMotherPi = RecoDecay::getMother(mcParticles, particle2, Pdg::kDStar, true, nullptr, 1); // last arguement 1 is written to consider immediate decay mother only
        auto indexMotherD0 = RecoDecay::getMother(mcParticles, particle1, Pdg::kDStar, true, nullptr, 1);
        bool correlationStatus = false;
        if (std::abs(particle2.pdgCode()) == kPiPlus && indexMotherPi >= 0 && indexMotherD0 >= 0 && indexMotherPi == indexMotherD0) {
          if (!storeAutoCorrelationFlag) {
            continue;
          }
          correlationStatus = true;
        }

        registry.fill(HIST("hTrackCounterGen"), 3); // fill after soft pion removal

        if (correlateD0WithLeadingParticle) {
          if (particle2.globalIndex() != leadingIndex) {
            continue;
          }
          registry.fill(HIST("hTrackCounterGen"), 4); // fill no. of tracks  have leading particle
        }

        auto getTracksSize = [&mcParticles](aod::McCollision const& /*collision*/) {
          int nTracks = 0;
          for (const auto& track : mcParticles) {
            if (track.isPhysicalPrimary() && std::abs(track.eta()) < 1.0) {
              nTracks++;
            }
          }
          return nTracks;
        };
        using BinningTypeMcGen = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::mccollision::PosZ, decltype(getTracksSize)>;
        BinningTypeMcGen corrBinningMcGen{{getTracksSize}, {zPoolBins, multPoolBinsMcGen}, true};
        int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), getTracksSize(mcCollision)));
        entryD0HadronPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                          particle2.eta() - particle1.eta(),
                          particle1.pt(),
                          particle2.pt(),
                          poolBin,
                          correlationStatus);
        entryD0HadronRecoInfo(massD0, massD0, 0); // dummy info
      }                                           // end inner loop (Tracks)
    }                                             // end outer loop (D0)
  }

  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processMcGen, "Process MC Gen mode", false);

  // ====================== Implement Event mixing on Data ===================================

  void processDataMixedEvent(SelectedCollisions const& collisions,
                             SelectedCandidatesData const& candidates,
                             SelectedTracks const& tracks)
  {
    for (const auto& collision : collisions) {
      registry.fill(HIST("hMultFT0M"), collision.multFT0M());
      registry.fill(HIST("hZvtx"), collision.posZ());
    }

    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelectedCollisions, SelectedCandidatesData, SelectedTracks, BinningType> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      // LOGF(info, "Mixed event collisions: Index = (%d, %d), tracks Size: (%d, %d), Z Vertex: (%f, %f), Pool Bin: (%d, %d)", c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size(), c1.posZ(), c2.posZ(), corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M())),corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()))); // For debug
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (std::abs(hfHelper.yD0(t1)) >= yCandMax) {
          continue;
        }

        // soft pion removal, signal status 1,3 for D0 and 2,3 for D0bar (SoftPi removed), signal status 11,13 for D0  and 12,13 for D0bar (only SoftPi)
        auto ePiK = RecoDecay::e(t1.pVectorProng0(), massPi) + RecoDecay::e(t1.pVectorProng1(), massK);
        auto eKPi = RecoDecay::e(t1.pVectorProng0(), massK) + RecoDecay::e(t1.pVectorProng1(), massPi);
        double invMassDstar1 = 0., invMassDstar2 = 0.;
        bool isSoftPiD0 = false, isSoftPiD0bar = false;
        auto pSum2 = RecoDecay::p2(t1.pVector(), t2.pVector());
        auto ePion = t2.energy(massPi);
        invMassDstar1 = std::sqrt((ePiK + ePion) * (ePiK + ePion) - pSum2);
        invMassDstar2 = std::sqrt((eKPi + ePion) * (eKPi + ePion) - pSum2);

        if (t1.isSelD0() >= selectionFlagD0) {
          if ((std::abs(invMassDstar1 - hfHelper.invMassD0ToPiK(t1)) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0 = true;
          }
        }
        if (t1.isSelD0bar() >= selectionFlagD0bar) {
          if ((std::abs(invMassDstar2 - hfHelper.invMassD0barToKPi(t1)) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0bar = true;
          }
        }

        int signalStatus = 0;
        if (t1.isSelD0() >= selectionFlagD0) {
          if (!isSoftPiD0) {
            signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0Only;
          } else {
            signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0OnlySoftPi;
          }
        }
        if (t1.isSelD0bar() >= selectionFlagD0bar) {
          if (!isSoftPiD0bar) {
            signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0barOnly;
          } else {
            signalStatus += aod::hf_correlation_d0_hadron::ParticleTypeData::D0barOnlySoftPi;
          }
        }
        bool correlationStatus = false;
        entryD0HadronPair(getDeltaPhi(t1.phi(), t2.phi()), t1.eta() - t2.eta(), t1.pt(), t2.pt(), poolBin, correlationStatus);
        entryD0HadronRecoInfo(hfHelper.invMassD0ToPiK(t1), hfHelper.invMassD0barToKPi(t1), signalStatus);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processDataMixedEvent, "Process data mixed event", false);

  // ====================== Implement Event mixing on McRec ===================================

  void processMcRecMixedEvent(SelectedCollisions const& collisions,
                              SelectedCandidatesMcRec const& candidates,
                              SelectedTracksMcRec const& tracks)
  {
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelectedCollisions, SelectedCandidatesMcRec, SelectedTracksMcRec, BinningType> pairMcRec{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    bool flagD0 = false;
    bool flagD0bar = false;
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcRec) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));

      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (std::abs(hfHelper.yD0(t1)) >= yCandMax) {
          continue;
        }

        // soft pion removal
        auto ePiK = RecoDecay::e(t1.pVectorProng0(), massPi) + RecoDecay::e(t1.pVectorProng1(), massK);
        auto eKPi = RecoDecay::e(t1.pVectorProng0(), massK) + RecoDecay::e(t1.pVectorProng1(), massPi);
        double invMassDstar1 = 0., invMassDstar2 = 0.;
        bool isSoftPiD0 = false, isSoftPiD0bar = false;
        auto pSum2 = RecoDecay::p2(t1.pVector(), t2.pVector());
        auto ePion = t2.energy(massPi);
        invMassDstar1 = std::sqrt((ePiK + ePion) * (ePiK + ePion) - pSum2);
        invMassDstar2 = std::sqrt((eKPi + ePion) * (eKPi + ePion) - pSum2);

        if (t1.isSelD0() >= selectionFlagD0) {
          if ((std::abs(invMassDstar1 - hfHelper.invMassD0ToPiK(t1)) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0 = true;
          }
        }

        if (t1.isSelD0bar() >= selectionFlagD0bar) {
          if ((std::abs(invMassDstar2 - hfHelper.invMassD0barToKPi(t1)) - softPiMass) < ptSoftPionMax) {
            isSoftPiD0bar = true;
          }
        }

        flagD0 = t1.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK);     // flagD0Signal 'true' if candidate1 matched to D0 (particle)
        flagD0bar = t1.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK); // flagD0Reflection 'true' if candidate1, selected as D0 (particle), is matched to D0bar (antiparticle)
        int signalStatus = 0;

        if (flagD0 && (t1.isSelD0() >= selectionFlagD0)) {
          if (!isSoftPiD0) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Sig); //  signalStatus += 1;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi); // signalStatus += 64;
          }
        } // signal case D0

        if (flagD0bar && (t1.isSelD0() >= selectionFlagD0)) {
          if (!isSoftPiD0) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Ref); //   signalStatus += 2;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi); // signalStatus += 64;
          }
        } // reflection case D0

        if (!flagD0 && !flagD0bar && (t1.isSelD0() >= selectionFlagD0)) {
          if (!isSoftPiD0) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0Bg); //  signalStatus += 4;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi);
          }
        } // background case D0

        if (flagD0bar && (t1.isSelD0bar() >= selectionFlagD0bar)) {
          if (!isSoftPiD0bar) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barSig); //  signalStatus += 8;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi);
          }
        } // signal case D0bar

        if (flagD0 && (t1.isSelD0bar() >= selectionFlagD0bar)) {
          if (!isSoftPiD0bar) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barRef); // signalStatus += 16;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi);
          }
        } // reflection case D0bar

        if (!flagD0 && !flagD0bar && (t1.isSelD0bar() >= selectionFlagD0bar)) {
          if (!isSoftPiD0bar) {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::D0barBg); //   signalStatus += 32;
          } else {
            SETBIT(signalStatus, aod::hf_correlation_d0_hadron::ParticleTypeMcRec::SoftPi);
          }
        } // background case D0bar
        registry.fill(HIST("hSignalStatusMERec"), signalStatus);
        bool correlationStatus = false;
        entryD0HadronPair(getDeltaPhi(t1.phi(), t2.phi()), t1.eta() - t2.eta(), t1.pt(), t2.pt(), poolBin, correlationStatus);
        entryD0HadronRecoInfo(hfHelper.invMassD0ToPiK(t1), hfHelper.invMassD0barToKPi(t1), signalStatus);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processMcRecMixedEvent, "Process Mixed Event MCRec", false);

  // ====================== Implement Event mixing on McGen ===================================

  void processMcGenMixedEvent(SelectedCollisionsMcGen const& collisions,
                              SelectedParticlesMcGen const& mcParticles)
  {

    auto getTracksSize = [&mcParticles, this](SelectedCollisionsMcGen::iterator const& collision) {
      int nTracks = 0;
      auto associatedTracks = mcParticles.sliceByCached(o2::aod::mcparticle::mcCollisionId, collision.globalIndex(), this->cache);
      for (const auto& track : associatedTracks) {
        if (track.isPhysicalPrimary() && std::abs(track.eta()) < 1.0) {
          nTracks++;
        }
      }
      return nTracks;
    };

    using BinningTypeMcGen = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::mccollision::PosZ, decltype(getTracksSize)>;
    BinningTypeMcGen corrBinningMcGen{{getTracksSize}, {zPoolBins, multPoolBinsMcGen}, true};

    auto tracksTuple = std::make_tuple(mcParticles, mcParticles);
    Pair<SelectedCollisionsMcGen, SelectedParticlesMcGen, SelectedParticlesMcGen, BinningTypeMcGen> pairMcGen{corrBinningMcGen, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairMcGen) {
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        // Check track t1 is D0
        if (std::abs(t1.pdgCode()) != Pdg::kD0) {
          continue;
        }

        double yD = RecoDecay::y(t1.pVector(), MassD0);
        if (std::abs(yD) >= yCandMax || t1.pt() <= ptCandMin || std::abs(t2.eta()) >= etaTrackMax || t2.pt() <= ptTrackMin) {
          continue;
        }
        if ((std::abs(t2.pdgCode()) != kElectron) && (std::abs(t2.pdgCode()) != kMuonMinus) && (std::abs(t2.pdgCode()) != kPiPlus) && (std::abs(t2.pdgCode()) != kKPlus) && (std::abs(t2.pdgCode()) != kProton)) {
          continue;
        }
        if (!t2.isPhysicalPrimary()) {
          continue;
        }

        // ==============================soft pion removal================================
        // method used: indexMother = -1 by default if the mother doesn't match with given PID of the mother. We find mother of pion if it is D* and mother of D0 if it is D*. If they are both positive and they both match each other, then it is detected as a soft pion

        auto indexMotherPi = RecoDecay::getMother(mcParticles, t2, Pdg::kDStar, true, nullptr, 1); // last arguement 1 is written to consider immediate decay mother only
        auto indexMotherD0 = RecoDecay::getMother(mcParticles, t1, Pdg::kDStar, true, nullptr, 1);
        if (std::abs(t2.pdgCode()) == kPiPlus && indexMotherPi >= 0 && indexMotherD0 >= 0 && indexMotherPi == indexMotherD0) {
          continue;
        }
        int poolBin = corrBinningMcGen.getBin(std::make_tuple(c2.posZ(), getTracksSize(c2)));
        bool correlationStatus = false;
        entryD0HadronPair(getDeltaPhi(t2.phi(), t1.phi()), t2.eta() - t1.eta(), t1.pt(), t2.pt(), poolBin, correlationStatus);
        entryD0HadronRecoInfo(massD0, massD0, 0); // dummy info
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorD0Hadrons, processMcGenMixedEvent, "Process Mixed Event MCGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCorrelatorD0HadronsSelection>(cfgc),
    adaptAnalysisTask<HfCorrelatorD0Hadrons>(cfgc)};
}
