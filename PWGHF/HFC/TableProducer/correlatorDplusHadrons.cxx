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

/// \file correlatorDplusHadrons.cxx
/// \author Shyam Kumar <shyam.kumar@cern.ch>

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
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies

double getDeltaPhi(double phiD, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PIHalf);
}

/// definition of variables for Dplus hadron pairs (in data-like, MC-reco and MC-kine tasks)
const int npTBinsMassAndEfficiency = o2::analysis::hf_cuts_dplus_to_pi_k_pi::nBinsPt;
std::vector<double> efficiencyDmeson(npTBinsMassAndEfficiency + 1);

// histogram binning definition
const int massAxisBins = 350;
const double massAxisMin = 1.7;
const double massAxisMax = 2.05;
const int phiAxisBins = 32;
const double phiAxisMin = -o2::constants::math::PIHalf;
const double phiAxisMax = 3. * o2::constants::math::PIHalf;
const int yAxisBins = 100;
const double yAxisMin = -2.;
const double yAxisMax = 2.;
const int ptDAxisBins = 180;
const double ptDAxisMin = 0.;
const double ptDAxisMax = 36.;

// definition of ME variables
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;

// Code to select a Dmeson in a collision
struct HfCorrelatorDplusHadronsDplusSelection {
  Produces<aod::DmesonSelection> dplusSel;

  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for Dplus"}; // 7 corresponds to topo+PID cuts
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  HfHelper hfHelper;
  SliceCache cache;

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>> selectedDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>> selectedDplusCandidatesMc = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  void processDplusSelectionData(aod::Collision const& collision,
                                 soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const& candidates)
  {
    bool isDplusFound = 0;
    if (selectedDplusCandidates.size() > 0) {
      auto selectedDplusCandidatesGrouped = selectedDplusCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

      for (const auto& candidate1 : selectedDplusCandidatesGrouped) {
        if (yCandMax >= 0. && std::abs(hfHelper.yDplus(candidate1)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
          continue;
        }
        isDplusFound = 1;
        break;
      }
    }
    dplusSel(isDplusFound);
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadronsDplusSelection, processDplusSelectionData, "Process Dplus Selection Data", false);

  void processDplusSelectionMcRec(aod::Collision const& collision,
                                  soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec> const& candidates)
  {
    bool isDplusFound = 0;
    if (selectedDplusCandidatesMc.size() > 0) {
      auto selectedDplusCandidatesMcGrouped = selectedDplusCandidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      for (const auto& candidate1 : selectedDplusCandidatesMcGrouped) {
        // check decay channel flag for candidate1
        if (!(candidate1.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(hfHelper.yDplus(candidate1)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
          continue;
        }
        isDplusFound = 1;
        break;
      }
    }
    dplusSel(isDplusFound);
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadronsDplusSelection, processDplusSelectionMcRec, "Process Dplus Selection MCRec", false);

  void processDplusSelectionMcGen(aod::McCollision const& mcCollision,
                                  aod::McParticles const& mcParticles)
  {
    bool isDplusFound = 0;
    for (const auto& particle1 : mcParticles) {
      if (std::abs(particle1.pdgCode()) != Pdg::kDPlus) {
        continue;
      }
      double yD = RecoDecay::y(std::array{particle1.px(), particle1.py(), particle1.pz()}, MassDPlus);
      if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle1.pt() < ptCandMin) {
        continue;
      }
      isDplusFound = 1;
      break;
    }
    dplusSel(isDplusFound);
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadronsDplusSelection, processDplusSelectionMcGen, "Process Dplus Selection MCGen", false);
};

/// Dplus-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfCorrelatorDplusHadrons {
  Produces<aod::DplusHadronPair> entryDplusHadronPair;
  Produces<aod::DplusHadronRecoInfo> entryDplusHadronRecoInfo;
  Produces<aod::Dplus> entryDplus;
  Produces<aod::Hadron> entryHadron;

  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for Dplus"}; // 7 corresponds to topo+PID cuts
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
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{efficiencyDmeson}, "Efficiency values for Dplus meson"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 100000.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis binsZVtx{"binsZVtx", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "Mixing bins - z-vertex"};
  ConfigurableAxis binsMultiplicityMc{"binsMultiplicityMc", {VARIABLE_WIDTH, 0.0f, 20.0f, 50.0f, 500.0f}, "Mixing bins - MC multiplicity"}; // In MCGen multiplicity is defined by counting tracks

  HfHelper hfHelper;
  SliceCache cache;
  BinningType corrBinning{{binsZVtx, binsMultiplicity}, true};

  // Event Mixing for the Data Mode
  using MySelCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::DmesonSelection>>;
  using MyTracks = soa::Filtered<aod::TracksWDca>;
  using MyCandidatesData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  // Event Mixing for the MCRec Mode
  using MyCandidatesMcRec = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  // Event Mixing for the MCGen Mode
  using McCollisionsSel = soa::Filtered<soa::Join<aod::McCollisions, aod::DmesonSelection>>;
  using McParticlesSel = soa::Filtered<aod::McParticles>;

  Filter collisionFilter = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (nabs(aod::track::pt) > ptTrackMin) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);
  Filter dplusFilter = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= 1;
  Filter collisionFilterGen = aod::hf_selection_dmeson_collision::dmesonSel == true;
  Filter particlesFilter = nabs(aod::mcparticle::pdgCode) == 411 || ((aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary);

  Preslice<aod::HfCand3Prong> perCol = aod::hf_cand::collisionId;

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>> selectedDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>> selectedDplusCandidatesMc = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "Dplus,Hadron candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0", "Dplus,Hadron candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1", "Dplus,Hadron candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng2", "Dplus,Hadron candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatus", "Dplus,Hadron candidates;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEta", "Dplus,Hadron candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhi", "Dplus,Hadron candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hY", "Dplus,Hadron candidates;candidate #it{#y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPtCandMCRec", "Dplus,Hadron candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0MCRec", "Dplus,Hadron candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1MCRec", "Dplus,Hadron candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng2MCRec", "Dplus,Hadron candidates - MC reco;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatusMCRec", "Dplus,Hadron candidates - MC reco;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEtaMCRec", "Dplus,Hadron candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCRec", "Dplus,Hadron candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCRec", "Dplus,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hMCEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMCGen", "Dplus,Hadron particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hEtaMCGen", "Dplus,Hadron particles - MC gen;particle #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCGen", "Dplus,Hadron particles - MC gen;particle #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCGen", "Dplus,Hadron candidates - MC gen;candidate #it{#y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hcountDplusHadronPerEvent", "Dplus,Hadron particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultV0M", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hZvtx", "z vertex;z vertex;entries", {HistType::kTH1F, {{200, -20., 20.}}}},
     {"hDplusBin", "Dplus selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}}},
     {"hTracksBin", "Tracks selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}}}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMassDplus_2D", "Dplus candidates;inv. mass (K^{-}#pi^{+}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDplusData", "Dplus candidates;inv. mass (K^{-}#pi^{+}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    registry.add("hMassDplusMCRec", "Dplus candidates;inv. mass (K^{-}#pi^{+}#pi^{+}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    registry.add("hMassDplusMCRecSig", "Dplus signal candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDplusMCRecBkg", "Dplus background candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hcountDplustriggersMCGen", "Dplus trigger particles - MC gen;;N of trigger Dplus", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    corrBinning = {{binsZVtx, binsMultiplicity}, true};
  }

  /// Dplus-hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision,
                   aod::TracksWDca const& tracks,
                   soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const& candidates, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int gCollisionId = collision.globalIndex();
    int64_t timeStamp = bc.timestamp();
    if (selectedDplusCandidates.size() > 0) {
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
          registry.fill(HIST("hTracksBin"), poolBin);
        }
      }

      registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
      if (nTracks < multMin || nTracks > multMax) {
        return;
      }
      registry.fill(HIST("hMultiplicity"), nTracks);

      auto selectedDplusCandidatesGrouped = selectedDplusCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      int cntDplus = 0;
      for (const auto& candidate1 : selectedDplusCandidatesGrouped) {
        if (yCandMax >= 0. && std::abs(hfHelper.yDplus(candidate1)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
          continue;
        }
        if (candidate1.pt() > ptTrackMax) {
          continue;
        }
        // check decay channel flag for candidate1
        if (!(candidate1.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi)) {
          continue;
        }
        double efficiencyWeight = 1.;
        if (applyEfficiency) {
          efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
        }
        // fill invariant mass plots and generic info from all Dplus candidates
        registry.fill(HIST("hMassDplus_2D"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        registry.fill(HIST("hMassDplusData"), hfHelper.invMassDplusToPiKPi(candidate1), efficiencyWeight);
        registry.fill(HIST("hPtCand"), candidate1.pt());
        registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
        registry.fill(HIST("hPtProng2"), candidate1.ptProng2());
        registry.fill(HIST("hEta"), candidate1.eta());
        registry.fill(HIST("hPhi"), RecoDecay::constrainAngle(candidate1.phi(), -o2::constants::math::PIHalf));
        registry.fill(HIST("hY"), hfHelper.yDplus(candidate1));
        registry.fill(HIST("hSelectionStatus"), candidate1.isSelDplusToPiKPi());
        registry.fill(HIST("hDplusBin"), poolBin);
        entryDplus(candidate1.phi(), candidate1.eta(), candidate1.pt(), hfHelper.invMassDplusToPiKPi(candidate1), poolBin, gCollisionId, timeStamp);
        // Dplus-Hadron correlation dedicated section
        // if the candidate is a Dplus, search for Hadrons and evaluate correlations
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          if (track.pt() < ptTrackMin) {
            continue;
          }
          if (std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
            continue; // Remove secondary tracks
          }
          // Removing Dplus daughters by checking track indices
          if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex()) || (candidate1.prong2Id() == track.globalIndex())) {
            continue;
          }
          entryDplusHadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                               track.eta() - candidate1.eta(),
                               candidate1.pt(),
                               track.pt(), poolBin);
          entryDplusHadronRecoInfo(hfHelper.invMassDplusToPiKPi(candidate1), 0);
          if (cntDplus == 0)
            entryHadron(track.phi(), track.eta(), track.pt(), poolBin, gCollisionId, timeStamp);
        } // Hadron Tracks loop
        cntDplus++;
      } // end outer Dplus loop
      registry.fill(HIST("hZvtx"), collision.posZ());
      registry.fill(HIST("hMultV0M"), collision.multFV0M());
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processData, "Process data", false);

  /// Dplus-Hadron correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRec(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision,
                    aod::TracksWDca const& tracks,
                    soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec> const& candidates)
  {
    if (selectedDplusCandidatesMc.size() > 0) {
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
          registry.fill(HIST("hTracksBin"), poolBin);
        }
      }
      registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
      if (nTracks < multMin || nTracks > multMax) {
        return;
      }
      registry.fill(HIST("hMultiplicity"), nTracks);

      auto selectedDplusCandidatesMcGrouped = selectedDplusCandidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      // MC reco level
      bool flagDplusSignal = false;
      for (const auto& candidate1 : selectedDplusCandidatesMcGrouped) {
        // check decay channel flag for candidate1
        if (!(candidate1.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(hfHelper.yDplus(candidate1)) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
          continue;
        }
        if (candidate1.pt() >= ptTrackMax) {
          continue;
        }
        double efficiencyWeight = 1.;
        if (applyEfficiency) {
          efficiencyWeight = 1. / efficiencyD->at(o2::analysis::findBin(binsPt, candidate1.pt()));
        }

        if (std::abs(candidate1.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi) {
          // fill per-candidate distributions from Dplus true candidates
          registry.fill(HIST("hPtCandMCRec"), candidate1.pt());
          registry.fill(HIST("hPtProng0MCRec"), candidate1.ptProng0());
          registry.fill(HIST("hPtProng1MCRec"), candidate1.ptProng1());
          registry.fill(HIST("hPtProng2MCRec"), candidate1.ptProng2());
          registry.fill(HIST("hEtaMCRec"), candidate1.eta());
          registry.fill(HIST("hPhiMCRec"), RecoDecay::constrainAngle(candidate1.phi(), -o2::constants::math::PIHalf));
          registry.fill(HIST("hYMCRec"), hfHelper.yDplus(candidate1));
          registry.fill(HIST("hSelectionStatusMCRec"), candidate1.isSelDplusToPiKPi());
        }
        // fill invariant mass plots from Dplus signal and background candidates
        registry.fill(HIST("hMassDplusMCRec"), hfHelper.invMassDplusToPiKPi(candidate1), efficiencyWeight);
        if (std::abs(candidate1.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi) { // also matched as Dplus
          registry.fill(HIST("hMassDplusMCRecSig"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassDplusMCRecBkg"), hfHelper.invMassDplusToPiKPi(candidate1), candidate1.pt(), efficiencyWeight);
        }
        registry.fill(HIST("hDplusBin"), poolBin);
        // Dplus-Hadron correlation dedicated section
        // if the candidate is selected as Dplus, search for Hadron and evaluate correlations
        flagDplusSignal = candidate1.flagMcMatchRec() == 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi;
        for (const auto& track : tracks) {
          if (std::abs(track.eta()) > etaTrackMax) {
            continue;
          }
          if (track.pt() < ptTrackMin) {
            continue;
          }
          // Removing Dplus daughters by checking track indices
          if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex()) || (candidate1.prong2Id() == track.globalIndex())) {
            continue;
          }
          if (std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
            continue; // Remove secondary tracks
          }
          entryDplusHadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                               track.eta() - candidate1.eta(),
                               candidate1.pt(),
                               track.pt(), poolBin);
          entryDplusHadronRecoInfo(hfHelper.invMassDplusToPiKPi(candidate1), flagDplusSignal);
        } // end inner loop (Tracks)

      } // end outer Dplus loop
      registry.fill(HIST("hZvtx"), collision.posZ());
      registry.fill(HIST("hMultV0M"), collision.multFV0M());
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processMcRec, "Process MC Reco mode", true);

  /// Dplus-Hadron correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::McCollision const& mcCollision,
                    aod::McParticles const& mcParticles)
  {
    int counterDplusHadron = 0;
    registry.fill(HIST("hMCEvtCount"), 0);

    auto getTracksSize = [&mcParticles](aod::McCollision const& collision) {
      int nTracks = 0;
      for (const auto& track : mcParticles) {
        if (track.isPhysicalPrimary() && std::abs(track.eta()) < 1.0) {
          nTracks++;
        }
      }
      return nTracks;
    };
    using BinningTypeMCGen = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::mccollision::PosZ, decltype(getTracksSize)>;
    BinningTypeMCGen corrBinningMcGen{{getTracksSize}, {binsZVtx, binsMultiplicityMc}, true};

    // MC gen level
    for (const auto& particle1 : mcParticles) {
      // check if the particle is Dplus  (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != Pdg::kDPlus) {
        continue;
      }
      double yD = RecoDecay::y(std::array{particle1.px(), particle1.py(), particle1.pz()}, MassDPlus);
      if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle1.pt() < ptCandMin) {
        continue;
      }
      registry.fill(HIST("hPtCandMCGen"), particle1.pt());
      registry.fill(HIST("hEtaMCGen"), particle1.eta());
      registry.fill(HIST("hPhiMCGen"), RecoDecay::constrainAngle(particle1.phi(), -o2::constants::math::PIHalf));
      registry.fill(HIST("hYMCGen"), yD);
      counterDplusHadron++;
      // Dplus Hadron correlation dedicated section
      // if it's a Dplus particle, search for Hadron and evaluate correlations
      if (std::abs(particle1.pdgCode()) != Pdg::kDPlus) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
        continue;
      }
      registry.fill(HIST("hcountDplustriggersMCGen"), 0, particle1.pt()); // to count trigger Dplus for normalisation)
      for (const auto& particle2 : mcParticles) {

        // Check Mother of particle 2
        bool flagMotherFound = false;
        for (const auto& m : particle2.mothers_as<aod::McParticles>()) {
          if (m.globalIndex() == particle1.globalIndex()) {
            flagMotherFound = true;
            break;
          }
        }
        if (flagMotherFound) {
          continue;
        }
        if (std::abs(particle2.eta()) > etaTrackMax) {
          continue;
        }
        if (particle2.pt() < ptTrackMin) {
          continue;
        }

        if ((std::abs(particle2.pdgCode()) != 11) && (std::abs(particle2.pdgCode()) != 13) && (std::abs(particle2.pdgCode()) != 211) && (std::abs(particle2.pdgCode()) != 321) && (std::abs(particle2.pdgCode()) != 2212)) {
          continue;
        }
        int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), getTracksSize(mcCollision)));
        entryDplusHadronPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                             particle2.eta() - particle1.eta(),
                             particle1.pt(),
                             particle2.pt(), poolBin);

      } // end inner loop
    }   // end outer loop
    registry.fill(HIST("hcountDplusHadronPerEvent"), counterDplusHadron);
    registry.fill(HIST("hZvtx"), mcCollision.posZ());
    registry.fill(HIST("hMultiplicity"), getTracksSize(mcCollision));
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processMcGen, "Process MC Gen mode", false);

  void processDataMixedEvent(MySelCollisions const& collisions,
                             MyCandidatesData const& candidates,
                             MyTracks const& tracks)
  {
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<MySelCollisions, MyCandidatesData, MyTracks, BinningType> pairData{corrBinning, 5, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      // LOGF(info, "Mixed event collisions: Index = (%d, %d), tracks Size: (%d, %d), Z Vertex: (%f, %f), Pool Bin: (%d, %d)", c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size(), c1.posZ(), c2.posZ(), corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFV0M())),corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFV0M()))); // For debug
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFV0M()));
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (yCandMax >= 0. && std::abs(hfHelper.yDplus(t1)) > yCandMax) {
          continue;
        }
        entryDplusHadronPair(getDeltaPhi(t1.phi(), t2.phi()), t1.eta() - t2.eta(), t1.pt(), t2.pt(), poolBin);
        entryDplusHadronRecoInfo(hfHelper.invMassDplusToPiKPi(t1), 0);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processDataMixedEvent, "Process Mixed Event Data", false);

  void processMcRecMixedEvent(MySelCollisions const& collisions,
                              MyCandidatesMcRec const& candidates,
                              MyTracks const& tracks)
  {
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<MySelCollisions, MyCandidatesMcRec, MyTracks, BinningType> pairMcRec{corrBinning, 5, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairMcRec) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFV0M()));
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (yCandMax >= 0. && std::abs(hfHelper.yDplus(t1)) > yCandMax) {
          continue;
        }
        entryDplusHadronPair(getDeltaPhi(t1.phi(), t2.phi()), t1.eta() - t2.eta(), t1.pt(), t2.pt(), poolBin);
        entryDplusHadronRecoInfo(hfHelper.invMassDplusToPiKPi(t1), 0);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processMcRecMixedEvent, "Process Mixed Event MCRec", false);

  void processMcGenMixedEvent(McCollisionsSel const& collisions,
                              McParticlesSel const& mcParticles)
  {

    auto getTracksSize = [&mcParticles, this](McCollisionsSel::iterator const& collision) {
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
    BinningTypeMcGen corrBinningMcGen{{getTracksSize}, {binsZVtx, binsMultiplicityMc}, true};

    auto tracksTuple = std::make_tuple(mcParticles, mcParticles);
    Pair<McCollisionsSel, McParticlesSel, McParticlesSel, BinningTypeMcGen> pairMcGen{corrBinningMcGen, 5, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairMcGen) {
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        // Check track t1 is Dplus
        if (std::abs(t1.pdgCode()) != Pdg::kDPlus) {
          continue;
        }

        double yD = RecoDecay::y(std::array{t1.px(), t1.py(), t1.pz()}, MassDPlus);
        if (yCandMax >= 0. && std::abs(yD) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && t1.pt() < ptCandMin) {
          continue;
        }

        if (std::abs(t2.eta()) > etaTrackMax) {
          continue;
        }
        if (t2.pt() < ptTrackMin) {
          continue;
        }
        int poolBin = corrBinningMcGen.getBin(std::make_tuple(c2.posZ(), getTracksSize(c2)));
        // LOGF(info, "Mixed event collisions: Index = (%d,%d), tracks Size: (%d,%d), Z Vertex: (%f), Pool Bin: (%d)", c1.globalIndex(), c2.globalIndex(), getTracksSize(c1), getTracksSize(c2), c2.posZ(), poolBin); // For debug
        entryDplusHadronPair(getDeltaPhi(t2.phi(), t1.phi()), t2.eta() - t1.eta(), t1.pt(), t2.pt(), poolBin);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorDplusHadrons, processMcGenMixedEvent, "Process Mixed Event MCGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDplusHadronsDplusSelection>(cfgc),
                      adaptAnalysisTask<HfCorrelatorDplusHadrons>(cfgc)};
}
