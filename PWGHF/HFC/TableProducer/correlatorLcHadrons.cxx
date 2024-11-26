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
/// Returns deltaPhi values in range [-pi/2., 3.*pi/2.], typically used for correlation studies
///
double getDeltaPhi(double phiLc, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiLc, -o2::constants::math::PIHalf);
}

/// definition of variables for Lc hadron pairs (in data-like, Mc-reco and Mc-kine tasks)
const int nBinsPtMassAndEfficiency = o2::analysis::hf_cuts_lc_to_p_k_pi::nBinsPt;
const double efficiencyLcDefault[nBinsPtMassAndEfficiency] = {};
auto vecEfficiencyLc = std::vector<double>{efficiencyLcDefault, efficiencyLcDefault + nBinsPtMassAndEfficiency};

// histogram binning definition
const int massAxisBins = 120;
const double massAxisMin = 1.98;
const double massAxisMax = 2.58;
const int phiAxisBins = 32;
const double phiAxisMin = -o2::constants::math::PIHalf;
const double phiAxisMax = 3. * o2::constants::math::PIHalf;
const int yAxisBins = 100;
const double yAxisMin = -2.;
const double yAxisMax = 2.;
const int ptLcAxisBins = 180;
const double ptLcAxisMin = 0.;
const double ptLcAxisMax = 36.;

// definition of ME variables
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;

using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::LcSelection>>;
using SelectedTracks = soa::Filtered<aod::TracksWDca>;
using SelectedCandidatesData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
using SelectedCandidatesMcRec = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>;
using SelectedCollisionsMcGen = soa::Filtered<soa::Join<aod::McCollisions, aod::LcSelection>>;
using SelectedTracksMcGen = soa::Filtered<aod::McParticles>;

// Code to select collisions with at least one Lambda_c
struct HfCorrelatorLcHadronsSelection {
  Produces<aod::LcSelection> lcSel;

  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  HfHelper hfHelper;
  SliceCache cache;

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelLc>> selectedLcCandidates = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>> selectedLcCandidatesMc = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc;

  // Returns false if the candidate does not pass cuts on decay type, y max, and pt min. Used for data and MC reco.
  template <typename T>
  bool kinematicCuts(const T& candidate)
  {
    // check decay channel flag for candidate
    if (!TESTBIT(candidate.hfflag(), aod::hf_cand_3prong::DecayType::LcToPKPi)) {
      return false;
    }
    if (yCandMax >= 0. && std::abs(hfHelper.yLc(candidate)) > yCandMax) {
      return false;
    }
    if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
      return false;
    }
    return true;
  }

  void processLcSelectionData(aod::Collision const& collision,
                              soa::Join<aod::HfCand3Prong, aod::HfSelLc> const&)
  {
    int isLcFound = 0;
    if (selectedLcCandidates.size() > 0) {
      auto selectedLcCandidatesGrouped = selectedLcCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

      for (const auto& candidate : selectedLcCandidatesGrouped) {
        if (!kinematicCuts(candidate)) {
          continue;
        }
        isLcFound = 1;
        break;
      }
    }
    lcSel(isLcFound);
  }
  PROCESS_SWITCH(HfCorrelatorLcHadronsSelection, processLcSelectionData, "Process Lc Collision Selection Data", true);

  void processLcSelectionMcRec(aod::Collision const& collision,
                               soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec> const&)
  {
    int isLcFound = 0;
    if (selectedLcCandidatesMc.size() > 0) {
      auto selectedLcCandidatesGroupedMc = selectedLcCandidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      for (const auto& candidate : selectedLcCandidatesGroupedMc) {
        if (!kinematicCuts(candidate)) {
          continue;
        }
        isLcFound = 1;
        break;
      }
    }
    lcSel(isLcFound);
  }
  PROCESS_SWITCH(HfCorrelatorLcHadronsSelection, processLcSelectionMcRec, "Process Lc Selection McRec", false);

  void processLcSelectionMcGen(aod::McCollision const&,
                               aod::McParticles const& mcParticles)
  {
    int isLcFound = 0;
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != Pdg::kLambdaCPlus) {
        continue;
      }
      double yL = RecoDecay::y(particle.pVector(), MassLambdaCPlus);
      if (yCandMax >= 0. && std::abs(yL) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
        continue;
      }
      isLcFound = 1;
      break;
    }
    lcSel(isLcFound);
  }
  PROCESS_SWITCH(HfCorrelatorLcHadronsSelection, processLcSelectionMcGen, "Process Lc Selection McGen", false);
};

// Lc-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via Mc truth)
struct HfCorrelatorLcHadrons {
  Produces<aod::LcHadronPair> entryLcHadronPair;
  Produces<aod::LcHadronRecoInfo> entryLcHadronRecoInfo;

  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying Lc efficiency weights"};
  Configurable<int> filterFlagLc{"filterFlagLc", true, "Flag for applying Collision Filter with Lc"};
  Configurable<int> filterFlagLcMc{"filterFlagLcMc", false, "Flag for applying Mc Collision Filter with Lc"};
  Configurable<int> nEventForMixedEvent{"nEventForMixedEvent", 5, "number of event to be mixed"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 0.0025, "max. DCAxy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 0.0025, "max. DCAz of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 50., "max. track pT"};
  Configurable<float> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<float> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyLc{"efficiencyLc", std::vector<double>{vecEfficiencyLc}, "Efficiency values for Lc"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 100000.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis binsZVtx{"binsZVtx", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "Mixing bins - z-vertex"};
  ConfigurableAxis binsMultiplicityMc{"binsMultiplicityMc", {VARIABLE_WIDTH, 0.0f, 20.0f, 50.0f, 500.0f}, "Mixing bins - MC multiplicity"}; // In MCGen multiplicity is defined by counting tracks
  Configurable<bool> storeAutoCorrelationFlag{"storeAutoCorrelationFlag", false, "Store flag that indicates if the track is paired to its D-meson mother instead of skipping it"};
  Configurable<bool> correlateLcWithLeadingParticle{"correlateLcWithLeadingParticle", false, "Switch for correlation of Lc baryons with leading particle only"};

  HfHelper hfHelper;
  SliceCache cache;
  BinningType corrBinning{{binsZVtx, binsMultiplicity}, true};
  int leadingIndex = 0;
  bool correlationStatus = false;

  // Filters for ME
  Filter collisionFilter = aod::hf_selection_lc_collision::lcSel >= filterFlagLc;
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (nabs(aod::track::pt) > ptTrackMin) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);
  Filter lcFilter = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= 1) || (aod::hf_sel_candidate_lc::isSelLcToPiKP >= 1);
  Filter collisionFilterGen = aod::hf_selection_lc_collision::lcSel >= filterFlagLcMc;
  Filter particlesFilter = nabs(aod::mcparticle::pdgCode) == 4122 || ((aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary);

  Preslice<aod::HfCand3Prong> perCol = aod::hf_cand::collisionId;

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelLc>> selectedLcCandidates = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>> selectedLcCandidatesMc = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc;

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "Lc,Hadron candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hPtProng0", "Lc,Hadron candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hPtProng1", "Lc,Hadron candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hPtProng2", "Lc,Hadron candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hSelectionStatusLcToPKPi", "Lc,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hSelectionStatusLcToPiKP", "Lc,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hEta", "Lc,Hadron candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhi", "Lc,Hadron candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hY", "Lc,Hadron candidates;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPtCandMcRec", "Lc,Hadron candidates - Mc reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hPtProng0McRec", "Lc,Hadron candidates - Mc reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hPtProng1McRec", "Lc,Hadron candidates - Mc reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hPtProng2McRec", "Lc,Hadron candidates - Mc reco;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hSelectionStatusMcRec", "Lc,Hadron candidates - Mc reco;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}}},
     {"hEtaMcRec", "Lc,Hadron candidates - Mc reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMcRec", "Lc,Hadron candidates - Mc reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMcRec", "Lc,Hadron candidates - Mc reco;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hMcEvtCount", "Event counter - Mc gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMcGen", "Lc,Hadron particles - Mc gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hPtParticleAssocMcRec", "Associated Particles - Mc Rec;Hadron #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hPtParticleAssocMcGen", "Associated Particles - Mc Gen;Hadron #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptLcAxisBins, ptLcAxisMin, ptLcAxisMax}}}},
     {"hEtaMcGen", "Lc,Hadron particles - Mc Gen;particle #it{#varphi};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMcGen", "Lc,Hadron particles - Mc Gen;particle #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMcGen", "Lc,Hadron candidates - Mc Gen;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hCountLcHadronPerEvent", "Lc,Hadron particles - Mc Gen;Number per event;entries", {HistType::kTH1F, {{21, -0.5, 20.5}}}},
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultT0M", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hZvtx", "z vertex;z vertex;entries", {HistType::kTH1F, {{200, -20., 20.}}}},
     {"hLcPoolBin", "Lc selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}}},
     {"hTracksPoolBin", "Tracks selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}}}}};

  void init(InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMassLcVsPt", "Lc candidates;inv. mass (p k #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassLcData", "Lc candidates;inv. mass (p k #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    registry.add("hMassLcMcRec", "Lc candidates;inv. mass (p k #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    registry.add("hMassLcVsPtMcRec", "Lc candidates;inv. mass (p k #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassLcMcRecSig", "Lc signal candidates - Mc reco;inv. mass (p k #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // mass histogram for Lc background candidates only
    registry.add("hMassLcMcRecBkg", "Lc background candidates - Mc reco;inv. mass (p k #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCountLctriggersMcGen", "Lc trigger particles - Mc gen;;N of trigger Lc", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    corrBinning = {{binsZVtx, binsMultiplicity}, true};
  }

  /// Lc-h correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via Mc truth)
  void processData(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision,
                   aod::TracksWDca const& tracks,
                   soa::Join<aod::HfCand3Prong, aod::HfSelLc> const&)
  {
    // protection against empty tables to be sliced
    if (selectedLcCandidates.size() == 0) {
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
        if (std::abs(track.eta()) > etaTrackMax) {
          continue;
        }
        if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
          continue;
        }
        nTracks++;
        registry.fill(HIST("hTracksPoolBin"), poolBin);
      }
    }
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < multMin || nTracks > multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    auto selectedLcCandidatesGrouped = selectedLcCandidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

    for (const auto& candidate : selectedLcCandidatesGrouped) {
      if (yCandMax >= 0. && std::abs(hfHelper.yLc(candidate)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
        continue;
      }
      if (candidate.pt() > ptTrackMax) {
        continue;
      }
      // check decay channel flag for candidate
      if (!TESTBIT(candidate.hfflag(), aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / efficiencyLc->at(o2::analysis::findBin(binsPt, candidate.pt()));
      }
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hPtProng2"), candidate.ptProng2());
      registry.fill(HIST("hEta"), candidate.eta());
      registry.fill(HIST("hPhi"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf));
      registry.fill(HIST("hY"), hfHelper.yLc(candidate));
      registry.fill(HIST("hLcPoolBin"), poolBin);
      if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
        registry.fill(HIST("hMassLcVsPt"), hfHelper.invMassLcToPKPi(candidate), candidate.pt(), efficiencyWeight);
        registry.fill(HIST("hMassLcData"), hfHelper.invMassLcToPKPi(candidate), efficiencyWeight);
        registry.fill(HIST("hSelectionStatusLcToPKPi"), candidate.isSelLcToPKPi());
      }
      if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
        registry.fill(HIST("hMassLcVsPt"), hfHelper.invMassLcToPiKP(candidate), candidate.pt(), efficiencyWeight);
        registry.fill(HIST("hMassLcData"), hfHelper.invMassLcToPiKP(candidate), efficiencyWeight);
        registry.fill(HIST("hSelectionStatusLcToPiKP"), candidate.isSelLcToPiKP());
      }
      // Lc-Hadron correlation dedicated section
      // if the candidate is a Lc, search for Hadrons and evaluate correlations

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

        // Remove Lc daughters by checking track indices
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
                            candidate.pt(),
                            track.pt(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPKPi(candidate), false);
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPiKP(candidate), false);
        }
      } // Hadron Tracks loop
    }   // end outer Lc loop
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultT0M"), collision.multFT0M());
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processData, "Process data", true);

  /// Lc-Hadron correlation process starts for McRec
  void processMcRec(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision,
                    aod::TracksWDca const& tracks,
                    soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec> const&)
  {
    if (selectedLcCandidatesMc.size() == 0) {
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
        if (std::abs(track.eta()) > etaTrackMax) {
          continue;
        }
        if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
          continue;
        }
        nTracks++;
        registry.fill(HIST("hTracksPoolBin"), poolBin);
      }
    }
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < multMin || nTracks > multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    auto selectedLcCandidatesGroupedMc = selectedLcCandidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

    // Mc reco level
    bool isLcSignal = false;
    for (const auto& candidate : selectedLcCandidatesGroupedMc) {
      // check decay channel flag for candidate
      if (!TESTBIT(candidate.hfflag(), aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(hfHelper.yLc(candidate)) > yCandMax) {
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
        efficiencyWeight = 1. / efficiencyLc->at(o2::analysis::findBin(binsPt, candidate.pt()));
      }
      isLcSignal = std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi;
      if (isLcSignal) {
        registry.fill(HIST("hPtCandMcRec"), candidate.pt());
        registry.fill(HIST("hPtProng0McRec"), candidate.ptProng0());
        registry.fill(HIST("hPtProng1McRec"), candidate.ptProng1());
        registry.fill(HIST("hPtProng2McRec"), candidate.ptProng2());
        registry.fill(HIST("hEtaMcRec"), candidate.eta());
        registry.fill(HIST("hPhiMcRec"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf));
        registry.fill(HIST("hYMcRec"), hfHelper.yLc(candidate));
        // LcToPKPi and LcToPiKP division
        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          registry.fill(HIST("hMassLcMcRec"), hfHelper.invMassLcToPKPi(candidate), efficiencyWeight);
          registry.fill(HIST("hMassLcMcRecSig"), hfHelper.invMassLcToPKPi(candidate), candidate.pt(), efficiencyWeight);
          registry.fill(HIST("hMassLcVsPtMcRec"), hfHelper.invMassLcToPKPi(candidate), candidate.pt(), efficiencyWeight);
          registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelLcToPKPi());
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          registry.fill(HIST("hMassLcMcRec"), hfHelper.invMassLcToPiKP(candidate), efficiencyWeight);
          registry.fill(HIST("hMassLcMcRecSig"), hfHelper.invMassLcToPiKP(candidate), candidate.pt(), efficiencyWeight);
          registry.fill(HIST("hMassLcVsPtMcRec"), hfHelper.invMassLcToPiKP(candidate), candidate.pt(), efficiencyWeight);
          registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelLcToPiKP());
        }
      } else {
        registry.fill(HIST("hPtCandMcRec"), candidate.pt());
        registry.fill(HIST("hPtProng0McRec"), candidate.ptProng0());
        registry.fill(HIST("hPtProng1McRec"), candidate.ptProng1());
        registry.fill(HIST("hPtProng2McRec"), candidate.ptProng2());
        registry.fill(HIST("hEtaMcRec"), candidate.eta());
        registry.fill(HIST("hPhiMcRec"), RecoDecay::constrainAngle(candidate.phi(), -o2::constants::math::PIHalf));
        registry.fill(HIST("hYMcRec"), hfHelper.yLc(candidate));
        // LcToPKPi and LcToPiKP division
        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          registry.fill(HIST("hMassLcMcRec"), hfHelper.invMassLcToPKPi(candidate), efficiencyWeight);
          registry.fill(HIST("hMassLcMcRecBkg"), hfHelper.invMassLcToPKPi(candidate), candidate.pt(), efficiencyWeight);
          registry.fill(HIST("hMassLcVsPtMcRec"), hfHelper.invMassLcToPKPi(candidate), candidate.pt(), efficiencyWeight);
          registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelLcToPKPi());
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          registry.fill(HIST("hMassLcMcRec"), hfHelper.invMassLcToPiKP(candidate), efficiencyWeight);
          registry.fill(HIST("hMassLcMcRecBkg"), hfHelper.invMassLcToPiKP(candidate), candidate.pt(), efficiencyWeight);
          registry.fill(HIST("hMassLcVsPtMcRec"), hfHelper.invMassLcToPiKP(candidate), candidate.pt(), efficiencyWeight);
          registry.fill(HIST("hSelectionStatusMcRec"), candidate.isSelLcToPiKP());
        }
      }
      registry.fill(HIST("hLcPoolBin"), poolBin);

      // Lc-Hadron correlation dedicated section
      // if the candidate is selected as Lc, search for Hadron ad evaluate correlations
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
        // Removing Lc daughters by checking track indices
        if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
          if (!storeAutoCorrelationFlag) {
            continue;
          }
          correlationStatus = true;
        }
        registry.fill(HIST("hPtParticleAssocMcRec"), track.pt());

        if (correlateLcWithLeadingParticle) {
          if (track.globalIndex() != leadingIndex) {
            continue;
          }
        }
        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPKPi(candidate), isLcSignal);
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                            track.eta() - candidate.eta(),
                            candidate.pt(),
                            track.pt(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPiKP(candidate), isLcSignal);
        }

      } // end inner loop (Tracks)
    }   // end outer Lc loop
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultT0M"), collision.multFT0M());
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processMcRec, "Process Mc Reco mode", false);

  /// Lc-Hadron correlation pair builder - for Mc gen-level analysis
  void processMcGen(aod::McCollision const& mcCollision,
                    soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticles)
  {
    int counterLcHadron = 0;
    registry.fill(HIST("hMcEvtCount"), 0);

    // find leading particle
    if (correlateLcWithLeadingParticle) {
      leadingIndex = findLeadingParticleMcGen(mcParticles, etaTrackMax.value, ptTrackMin.value);
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
    BinningTypeMcGen corrBinningMcGen{{getTracksSize}, {binsZVtx, binsMultiplicityMc}, true};

    // Mc gen level
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != Pdg::kLambdaCPlus) {
        continue;
      }
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {
        double yL = RecoDecay::y(particle.pVector(), MassLambdaCPlus);
        if (yCandMax >= 0. && std::abs(yL) > yCandMax) {
          continue;
        }
        if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
          continue;
        }
        registry.fill(HIST("hPtCandMcGen"), particle.pt());
        registry.fill(HIST("hEtaMcGen"), particle.eta());
        registry.fill(HIST("hPhiMcGen"), RecoDecay::constrainAngle(particle.phi(), -o2::constants::math::PIHalf));
        registry.fill(HIST("hYMcGen"), yL);
        counterLcHadron++;

        for (const auto& particleAssoc : mcParticles) {
          bool flagMotherFound = false;
          for (const auto& m : particleAssoc.mothers_as<aod::McParticles>()) {
            if (m.globalIndex() == particle.globalIndex()) {
              flagMotherFound = true;
              break;
            }
          }
          if (flagMotherFound) {
            continue;
          }
          if (std::abs(particleAssoc.eta()) > etaTrackMax) {
            continue;
          }
          if (particleAssoc.pt() < ptTrackMin) {
            continue;
          }

          if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particle.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
            if (!storeAutoCorrelationFlag) {
              continue;
            }
            correlationStatus = true;
          }

          if (correlateLcWithLeadingParticle) {
            if (particleAssoc.globalIndex() != leadingIndex) {
              continue;
            }
          }

          int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), getTracksSize(mcCollision)));
          registry.fill(HIST("hPtParticleAssocMcGen"), particleAssoc.pt());
          entryLcHadronPair(getDeltaPhi(particleAssoc.phi(), particle.phi()),
                            particleAssoc.eta() - particle.eta(),
                            particle.pt(),
                            particleAssoc.pt(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(MassLambdaCPlus, true);
        } // end inner loop
      }
    } // end outer loop
    registry.fill(HIST("hCountLcHadronPerEvent"), counterLcHadron);
    registry.fill(HIST("hZvtx"), mcCollision.posZ());
    registry.fill(HIST("hMultiplicity"), getTracksSize(mcCollision));
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processMcGen, "Process Mc Gen mode", false);

  void processDataMixedEvent(SelectedCollisions const& collisions,
                             SelectedCandidatesData const& candidates,
                             SelectedTracks const& tracks)
  {
    if (candidates.size() == 0) {
      return;
    }
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelectedCollisions, SelectedCandidatesData, SelectedTracks, BinningType> pairData{corrBinning, nEventForMixedEvent, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!TESTBIT(t1.hfflag(), aod::hf_cand_3prong::DecayType::LcToPKPi)) {
          continue;
        }
        if (yCandMax >= 0. && std::abs(hfHelper.yLc(t1)) > yCandMax) {
          continue;
        }
        // LcToPKPi and LcToPiKP division
        if (t1.isSelLcToPKPi() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(t1.phi(), t2.phi()),
                            t1.eta() - t2.eta(),
                            t1.pt(),
                            t2.pt(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPKPi(t1), false);
        }
        if (t1.isSelLcToPiKP() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(t1.phi(), t2.phi()),
                            t1.eta() - t2.eta(),
                            t1.pt(),
                            t2.pt(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPiKP(t1), false);
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processDataMixedEvent, "Process Mixed Event Data", false);

  void processMcRecMixedEvent(SelectedCollisions const& collisions,
                              SelectedCandidatesMcRec const& candidates,
                              SelectedTracks const& tracks)
  {
    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<SelectedCollisions, SelectedCandidatesMcRec, SelectedTracks, BinningType> pairMcRec{corrBinning, nEventForMixedEvent, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairMcRec) {
      int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (yCandMax >= 0. && std::abs(hfHelper.yLc(t1)) > yCandMax) {
          continue;
        }
        if (t1.isSelLcToPKPi() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(t1.phi(), t2.phi()),
                            t1.eta() - t2.eta(),
                            t1.pt(),
                            t2.pt(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPKPi(t1), false);
        }
        if (t1.isSelLcToPiKP() >= selectionFlagLc) {
          entryLcHadronPair(getDeltaPhi(t1.phi(), t2.phi()),
                            t1.eta() - t2.eta(),
                            t1.pt(),
                            t2.pt(),
                            poolBin,
                            correlationStatus);
          entryLcHadronRecoInfo(hfHelper.invMassLcToPiKP(t1), false);
        }
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorLcHadrons, processMcRecMixedEvent, "Process Mixed Event McRec", false);

  void processMcGenMixedEvent(SelectedCollisionsMcGen const& collisions,
                              SelectedTracksMcGen const& mcParticles)
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
    BinningTypeMcGen corrBinningMcGen{{getTracksSize}, {binsZVtx, binsMultiplicityMc}, true};

    auto tracksTuple = std::make_tuple(mcParticles, mcParticles);
    Pair<SelectedCollisionsMcGen, SelectedTracksMcGen, SelectedTracksMcGen, BinningTypeMcGen> pairMcGen{corrBinningMcGen, nEventForMixedEvent, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairMcGen) {
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        // Check track t1 is Lc
        if (std::abs(t1.pdgCode()) != Pdg::kLambdaCPlus) {
          continue;
        }

        double yL = RecoDecay::y(t1.pVector(), MassLambdaCPlus);
        if (yCandMax >= 0. && std::abs(yL) > yCandMax) {
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
        entryLcHadronPair(getDeltaPhi(t1.phi(), t2.phi()),
                          t1.eta() - t2.eta(),
                          t1.pt(),
                          t2.pt(),
                          poolBin,
                          correlationStatus);
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
