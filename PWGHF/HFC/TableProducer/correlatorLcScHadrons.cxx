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

/// \file correlatorLcScHadrons.cxx
/// \brief Lc-Hadrons correlator task - data-like, Mc-Reco and Mc-Gen analyses
///
/// \author Marianna Mazzilli <marianna.mazzilli@cern.ch>
/// \author Zhen Zhang <zhenz@cern.ch>
/// \author Ravindra Singh <ravindra.singh@cern.ch>

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFC/Utils/utilsCorrelations.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TPDGCode.h>
#include <TRandom3.h>

#include <Rtypes.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>

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
struct HfCorrelatorLcScHadronsSelection {
  Produces<aod::LcSelection> candSel;

  Configurable<bool> useSel8{"useSel8", true, "Flag for applying sel8 for collision selection"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<bool> doSelLcCollision{"doSelLcCollision", true, "Select collisions with at least one Lc"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  HfHelper hfHelper;
  SliceCache cache;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CandsLcDataFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using CandsLcMcRecFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>;
  using CandsScMcRec = soa::Join<aod::HfCandSc, aod::HfCandScMcRec>;
  using CandidatesLcMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;
  using CandidatesScMcGen = soa::Join<aod::McParticles, aod::HfCandScMcGen>;
  // filter on selection of Lc and decay channel Lc->PKPi
  Filter lcFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);

  template <bool isCandSc, typename CollType, typename CandType>
  void selectionCollision(CollType const& collision, CandType const& candidates)
  {
    bool isSelColl = true;
    bool isCandFound = false;
    bool isSel8 = true;
    bool isNosameBunchPileUp = true;
    double yCand = -999.;
    const int chargeScZero = 0;
    if (doSelLcCollision) {
      for (const auto& candidate : candidates) {

        if constexpr (isCandSc) {
          int8_t chargeCand = candidate.charge();

          if (chargeCand == chargeScZero) {
            yCand = hfHelper.ySc0(candidate);
          } else {
            yCand = hfHelper.yScPlusPlus(candidate);
          }

        } else {
          yCand = hfHelper.yLc(candidate);
        }

        if (std::abs(yCand) > yCandMax || candidate.pt() < ptCandMin) {
          isCandFound = false;
          continue;
        }
        isCandFound = true;
        break;
      }
    }
    if (useSel8) {
      isSel8 = collision.sel8();
    }
    if (selNoSameBunchPileUpColl) {
      isNosameBunchPileUp = static_cast<bool>(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup));
    }
    isSelColl = isCandFound && isSel8 && isNosameBunchPileUp;
    candSel(isSelColl);
  }

  template <bool isCandSc, typename CandType>
  void selectionCollisionMcGen(CandType const& mcParticles)
  {
    bool isCandFound = false;
    double massCand = -999.0;
    for (const auto& particle : mcParticles) {

      isCandFound = matchCandAndMass<isCandSc>(particle, massCand);
      if (!isCandFound) {
        continue;
      }

      double yCand = RecoDecay::y(particle.pVector(), massCand);
      if (std::abs(yCand) > yCandMax || particle.pt() < ptCandMin) {
        isCandFound = false;
        continue;
      }

      isCandFound = true;
      break;
    }
    candSel(isCandFound);
  }

  /// Code to select collisions with at least one Lc - for real data and data-like analysis
  void processLcSelection(SelCollisions::iterator const& collision,
                          CandsLcDataFiltered const& candidates)
  {
    selectionCollision<false>(collision, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadronsSelection, processLcSelection, "Process Lc Collision Selection for Data and Mc", true);

  void processScSelection(SelCollisions::iterator const& collision,
                          aod::HfCandSc const& candidates)
  {
    selectionCollision<true>(collision, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadronsSelection, processScSelection, "Process Sc Collision Selection for Data and Mc", false);

  void processLcSelectionMcRec(SelCollisions::iterator const& collision,
                               CandsLcMcRecFiltered const& candidates)
  {
    selectionCollision<false>(collision, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadronsSelection, processLcSelectionMcRec, "Process Lc Selection McRec", false);

  void processScSelectionMcRec(SelCollisions::iterator const& collision,
                               CandsScMcRec const& candidates)
  {
    selectionCollision<true>(collision, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadronsSelection, processScSelectionMcRec, "Process Sc Selection McRec", false);

  void processLcSelectionMcGen(aod::McCollision const&,
                               CandidatesLcMcGen const& mcParticles)
  {
    selectionCollisionMcGen<false>(mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadronsSelection, processLcSelectionMcGen, "Process Lc Selection McGen", false);

  void processScSelectionMcGen(aod::McCollision const&,
                               CandidatesScMcGen const& mcParticles)
  {
    selectionCollisionMcGen<true>(mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadronsSelection, processScSelectionMcGen, "Process Lc Selection McGen", false);
};

// Lc-Hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via Mc truth)
struct HfCorrelatorLcScHadrons {
  Produces<aod::LcHadronPair> entryCandHadronPair;
  Produces<aod::LcHadronPairY> entryCandHadronPairY;
  Produces<aod::LcHadronPairTrkPID> entryCandHadronPairTrkPID;
  Produces<aod::LcHadronRecoInfo> entryCandHadronRecoInfo;
  Produces<aod::LcHadronMlInfo> entryCandHadronMlInfo;
  Produces<aod::LcRecoInfo> entryCandCandRecoInfo;
  Produces<aod::LcHadronGenInfo> entryCandHadronGenInfo;
  Produces<aod::LcGenInfo> entryCandCandGenInfo;
  Produces<aod::TrkRecInfoLc> entryTrackRecoInfo;
  Produces<aod::Lc> entryCand;
  Produces<aod::Hadron> entryHadron;
  Produces<aod::LcHadronTrkPID> entryTrkPID;
  Produces<aod::CandChargePair> entryPairCandCharge;
  Produces<aod::CandCharge> entryCandCharge;

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
  Configurable<bool> calEffEventWithCand{"calEffEventWithCand", true, "Calculate the efficiency of Lc candidate"};
  Configurable<float> eventFractionToAnalyze{"eventFractionToAnalyze", -1, "Fraction of events to analyze (use only for ME offline on very large samples)"};

  HfHelper hfHelper;
  SliceCache cache;
  Service<o2::framework::O2DatabasePDG> pdg;
  int leadingIndex = 0;
  int poolBin = 0;
  int poolBinLc = 0;
  bool correlationStatus = false;
  bool isPrompt = false;
  bool isNonPrompt = false;
  bool isSignal = false;

  TRandom3* rnd = new TRandom3(0);
  // std::vector<float> outputMl = {-1., -1., -1.};
  std::vector<float> outputMlPKPi = {-1., -1., -1.};
  std::vector<float> outputMlPiKP = {-1., -1., -1.};

  // Event Mixing for the Data Mode
  // using SelCollisionsWithSc = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
  using SelCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::LcSelection>>;
  using SelCollisionsMc = soa::Filtered<soa::Join<aod::McCollisions, aod::LcSelection, aod::MultsExtraMC>>; // collisionFilter applied

  using CandsLcData = soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>;
  using CandsLcDataFiltered = soa::Filtered<CandsLcData>;

  // Event Mixing for the MCRec Mode
  using CandsLcMcRec = soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc, aod::HfMlLcToPKPi>;
  using CandsLcMcRecFiltered = soa::Filtered<CandsLcMcRec>;
  using CandidatesLcMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>; // flagLcFilter applied
  using CandsScMcRec = soa::Join<aod::HfCandSc, aod::HfCandScMcRec>;
  using CandidatesScMcGen = soa::Join<aod::McParticles, aod::HfCandScMcGen>;
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

  Preslice<aod::McParticles> perTrueCollision = o2::aod::mcparticle::mcCollisionId;
  Preslice<aod::Tracks> perCollisionID = aod::track::collisionId;
  Preslice<aod::HfCand3Prong> cand3ProngPerCol = aod::hf_cand::collisionId;
  Preslice<aod::HfCandSc> csndScPerCol = aod::hf_cand::collisionId;

  // configurable axis definition
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 100000.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis binsZVtx{"binsZVtx", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "Mixing bins - z-vertex"};
  ConfigurableAxis binsMultiplicityMc{"binsMultiplicityMc", {VARIABLE_WIDTH, 0.0f, 20.0f, 50.0f, 500.0f}, "Mixing bins - MC multiplicity"}; // In MCGen multiplicity is defined by counting tracks
  ConfigurableAxis binsBdtScore{"binsBdtScore", {100, 0., 1.}, "Bdt output scores"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 6000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsCandMass{"binsCandMass", {200, 1.98, 2.58}, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};

  BinningType corrBinning{{binsZVtx, binsMultiplicity}, true};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext&)
  {
    AxisSpec axisCandMass = {binsCandMass, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};
    AxisSpec axisEta = {binsEta, "#it{eta}"};
    AxisSpec axisPhi = {binsPhi, "#it{#varphi}"};
    AxisSpec axisPtLc = {static_cast<std::vector<double>>(binsPtLc), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {static_cast<std::vector<double>>(binsPtHadron), "#it{p}_{T} Hadron (GeV/#it{c})"};
    AxisSpec axisPtTrack = {500, 0, 50, "#it{p}_{T} Hadron (GeV/#it{c})"};
    AxisSpec axisMultiplicity = {binsMultiplicity, "Multiplicity"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisPosZ = {binsZVtx, "PosZ"};
    AxisSpec axisBdtScore = {binsBdtScore, "Bdt score"};
    AxisSpec axisPoolBin = {binsPoolBin, "PoolBin"};
    AxisSpec axisRapidity = {100, -2, 2, "Rapidity"};
    AxisSpec axisSign = {5, -2.5, 2.5, "Sign"};

    registry.add("hPtCand", "Lc,Hadron candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtProng0", "Lc,Hadron candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtProng1", "Lc,Hadron candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtProng2", "Lc,Hadron candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hSelectionStatusLcToPKPi", "Lc,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}});
    registry.add("hSelectionStatusLcToPiKP", "Lc,Hadron candidates;selection status;entries", {HistType::kTH1F, {{8, -0.5, 7.5}}});
    registry.add("hEta", "Lc,Hadron candidates;candidate #it{#eta};entries", {HistType::kTH1F, {axisEta}});
    registry.add("hPhi", "Lc,Hadron candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {axisPhi}});
    registry.add("hcountCandHadronPerEvent", "Lc,Hadron particles - MC gen;Number per event;entries", {HistType::kTH1F, {{21, -0.5, 20.5}}});
    registry.add("hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hMultFT0M", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hZvtx", "z vertex;z vertex;entries", {HistType::kTH1F, {{200, -20., 20.}}});
    registry.add("hCandBin", "Lc selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}});
    registry.add("hTracksBin", "Tracks selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}});
    registry.add("hMassLcVsPt", "Lc candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisCandMass}, {axisPtLc}}});
    registry.add("hMassScVsPtVsSign", "Sc candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});sign;entries", {HistType::kTH3F, {{axisCandMass}, {axisPtLc}, {axisSign}}});
    registry.add("hMassLcData", "Lc candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{axisCandMass}}});
    registry.add("hLcPoolBin", "Lc candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hTracksPoolBin", "Particles associated pool bin", {HistType::kTH1F, {axisPoolBin}});
    // Histograms for MC Reco analysis
    registry.add("hMcEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}});
    registry.add("hMassLcMcRecBkg", "Lc background candidates - MC reco;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisCandMass}, {axisPtLc}}});
    registry.add("hPtCandSig", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtCandSigPrompt", "Lc,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtCandSigNonPrompt", "Lc,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtCandMcRecBkg", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtLc}});
    registry.add("hEtaSig", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiSig", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
    registry.add("hY", "Lc,Hadron candidates;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hYSig", "Lc,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hPtCandMcRecSigPrompt", "Lc,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtCandMcRecSigNonPrompt", "Lc,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtLc}});
    registry.add("hEtaMcRecBkg", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcRecBkg", "Lc,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
    registry.add("hYMcRecBkg", "Lc,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hFakeTracksMcRec", "Fake tracks - MC Rec", {HistType::kTH1F, {axisPtHadron}});
    registry.add("hPtParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtLc}}});
    registry.add("hPtTracksVsSignRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtTrack}, {axisSign}}});
    registry.add("hPtTracksVsSignRecTrue", "Associated Particle - MC Rec (True)", {HistType::kTH2F, {{axisPtTrack}, {axisSign}}});
    registry.add("hPtTracksVsSignGen", "Associated Particle - MC Gen", {HistType::kTH2F, {{axisPtTrack}, {axisSign}}});
    registry.add("hPtPrimaryParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtLc}}});
    registry.add("hPtVsMultiplicityMcRecPrompt", "Multiplicity FT0M - MC Rec Prompt", {HistType::kTH2F, {{axisPtLc}, {axisMultFT0M}}});
    registry.add("hPtVsMultiplicityMcRecNonPrompt", "Multiplicity FT0M - MC Rec Non Prompt", {HistType::kTH2F, {{axisPtLc}, {axisMultFT0M}}});
    // Histograms for MC Gen analysis
    registry.add("hcountCandtriggersMcGen", "Lc trigger particles - MC gen;;N of trigger Lc", {HistType::kTH2F, {{1, -0.5, 0.5}, {axisPtLc}}});
    registry.add("hPtCandMcGen", "Lc,Hadron particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtLc}});
    registry.add("hYMcGen", "Lc,Hadron candidates - MC gen;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hPtCandMcGenPrompt", "Lc,Hadron particles - MC Gen Prompt", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtCandVsChargeMcGenPrompt", "Charm Hadron particles - MC Gen Prompt", {HistType::kTH2F, {{axisPtLc}, {axisSign}}});
    registry.add("hPtCandMcGenNonPrompt", "Charm Hadron particles - MC Gen Non Prompt", {HistType::kTH1F, {axisPtLc}});
    registry.add("hPtCandVsChargeMcGenNonPrompt", "Lc,Hadron particles - MC Gen Non Prompt", {HistType::kTH2F, {{axisPtLc}, {axisSign}}});
    registry.add("hPtParticleAssocMcGen", "Associated Particle - MC Gen", {HistType::kTH1F, {axisPtHadron}});
    registry.add("hEtaMcGen", "Lc,Hadron particles - MC Gen", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcGen", "Lc,Hadron particles - MC Gen", {HistType::kTH1F, {axisPhi}});
    registry.add("hMultFT0AMcGen", "Lc,Hadron multiplicity FT0A - MC Gen", {HistType::kTH1F, {axisMultiplicity}});

    corrBinning = {{binsZVtx, binsMultiplicity}, true};
  }

  template <typename MlProbType>
  void fillMlOutput(MlProbType const& mlProb, std::vector<float>& outputMl)
  {
    for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
      outputMl[iclass] = mlProb[classMl->at(iclass)];
    }
  };

  template <bool isCandSc, typename CandType>
  double estimateY(CandType const& candidate)
  {
    double y = -999.;
    const int chargeScZero = 0;
    if constexpr (isCandSc) {
      int8_t chargeCand = candidate.charge();

      if (chargeCand == chargeScZero) {
        y = hfHelper.ySc0(candidate);
      } else {
        y = hfHelper.yScPlusPlus(candidate);
      }

    } else {
      y = hfHelper.yLc(candidate);
    }
    return y;
  }

  template <typename T1, typename T2, typename McPart>
  void calculateTrkEff(T1 const& trackPos1, T2 const& trackPos2, McPart const& mcParticles)
  {
    // genrated tracks
    decltype(trackPos1.template mcParticle_as<aod::McParticles>()) mctrk{};
    if (trackPos1.has_mcParticle()) { // ambiguous tracks should be small
      mctrk = trackPos1.template mcParticle_as<aod::McParticles>();
    } else if (trackPos2.has_mcParticle()) {
      mctrk = trackPos2.template mcParticle_as<aod::McParticles>();
    } else {
      return;
    }
    auto gentracks = mcParticles.sliceBy(perTrueCollision, mctrk.mcCollisionId());
    for (const auto& track : gentracks) {
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

      auto motherTrkGen = mcParticles.iteratorAt(track.mothersIds()[0]);
      if (std::abs(motherTrkGen.pdgCode()) == kLambdaCPlus)
        continue;

      auto chargeTrack = pdg->GetParticle(track.pdgCode())->Charge(); // Retrieve charge
      registry.fill(HIST("hPtTracksVsSignGen"), track.pt(), chargeTrack / (std::abs(chargeTrack)));
    }
  }
  template <bool isMcRec, typename TrackType, typename CandType, typename McPart>
  void fillCorrelationTable(bool trkPidFill, TrackType const& track, CandType const& candidate,
                            const std::vector<float>& outMl, int binPool, int8_t correlStatus,
                            double yCand, int signCand, McPart const& mcParticles)
  {
    bool isPhysicalPrimary = false;
    int trackOrigin = -1;

    entryCandHadronPair(getDeltaPhi(track.phi(), candidate.phi()),
                        track.eta() - candidate.eta(),
                        candidate.pt(),
                        track.pt() * track.sign(),
                        binPool,
                        correlStatus);
    entryCandHadronPairY(track.rapidity(MassProton) - yCand);
    entryCandHadronMlInfo(outMl[0], outMl[1]);
    entryTrackRecoInfo(track.dcaXY(), track.dcaZ(), track.tpcNClsCrossedRows());
    entryPairCandCharge(signCand);
    if (trkPidFill) {
      entryCandHadronPairTrkPID(track.tpcNSigmaPr(), track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tofNSigmaPr(), track.tofNSigmaKa(), track.tofNSigmaPi());
    }
    if constexpr (isMcRec) {
      if (track.has_mcParticle()) {
        auto mcParticle = track.template mcParticle_as<aod::McParticles>();
        isPhysicalPrimary = mcParticle.isPhysicalPrimary();
        trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
        entryCandHadronGenInfo(isPrompt, isPhysicalPrimary, trackOrigin);
      } else {
        entryCandHadronGenInfo(isPrompt, false, 0);
        registry.fill(HIST("hFakeTracksMcRec"), track.pt());
      }
      registry.fill(HIST("hPtParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
      if (isPhysicalPrimary) {
        registry.fill(HIST("hPtPrimaryParticleAssocVsCandMcRec"), track.pt(), candidate.pt());
      }
    }
  }

  template <bool isMcRec, bool isCandSc, typename CollisionType, typename CandType, typename TrackType>
  void doSameEvent(CollisionType const& collision,
                   TrackType const& tracks,
                   CandType const& candidates,
                   aod::McParticles const* mcParticles = nullptr)
  {

    int nTracks = 0;
    int64_t timeStamp = 0;
    const int chargeScPlusPlus = 2;
    const int chargeScZero = 0;
    bool skipMixedEventTableFilling = false;
    float multiplicityFT0M = collision.multFT0M();
    int gCollisionId = collision.globalIndex();
    if (candidates.size() == 0) {
      return;
    }

    if (eventFractionToAnalyze > 0) {
      if (rnd->Uniform(0, 1) > eventFractionToAnalyze) {
        skipMixedEventTableFilling = true;
      }
    }

    if constexpr (!isMcRec) {
      timeStamp = collision.template bc_as<aod::BCsWithTimestamps>().timestamp();
    }

    poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), multiplicityFT0M));
    if (correlateLcWithLeadingParticle) {
      leadingIndex = findLeadingParticle(tracks, etaTrackMax.value);
    }

    // Count good tracks
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

    int countCand = 1;

    for (const auto& candidate : candidates) {

      int8_t chargeCand = 3;
      double efficiencyWeightCand = 1.;
      double yCand = -999.0;
      double etaCand = -999.0;
      double ptCandLc = -999.0;
      double ptCand = -999.0;
      double phiCand = -999.0;
      double massCandPKPi = -999.0;
      double massCandPiKP = -999.0;
      bool selLcPKPi = false;
      bool selLcPiKP = false;

      yCand = estimateY<isCandSc>(candidate);
      etaCand = candidate.eta();
      ptCand = candidate.pt();
      phiCand = RecoDecay::constrainAngle(candidate.phi(), -PIHalf);

      if ((std::abs(yCand) > yCandMax) || ptCand < ptCandMin || ptCand > ptCandMax) {
        continue;
      }

      registry.fill(HIST("hY"), yCand);
      registry.fill(HIST("hPtCand"), ptCand);
      registry.fill(HIST("hEta"), etaCand);
      registry.fill(HIST("hPhi"), phiCand);
      registry.fill(HIST("hCandBin"), poolBin);

      if (applyEfficiency) {
        efficiencyWeightCand = 1. / efficiencyLc->at(o2::analysis::findBin(binsPtEfficiencyLc, ptCand));
      }

      if constexpr (isMcRec) {
        isPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
        isNonPrompt = candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
      }

      if constexpr (isCandSc) {
        chargeCand = candidate.charge();
        const auto& candidateLc = candidate.template prongLc_as<CandsLcData>();
        ptCandLc = candidateLc.pt();
        selLcPKPi = (candidateLc.isSelLcToPKPi() >= selectionFlagLc) && (candidate.statusSpreadLcMinvPKPiFromPDG());
        selLcPiKP = (candidateLc.isSelLcToPiKP() >= selectionFlagLc) && (candidate.statusSpreadLcMinvPiKPFromPDG());
        if (selLcPKPi) {
          const auto& probs = candidateLc.mlProbLcToPKPi();
          fillMlOutput(probs, outputMlPKPi);
          massCandPKPi = hfHelper.invMassScRecoLcToPKPi(candidate, candidateLc);
        }
        if (selLcPiKP) {
          const auto& probs = candidateLc.mlProbLcToPiKP();
          fillMlOutput(probs, outputMlPiKP);
          massCandPiKP = hfHelper.invMassScRecoLcToPiKP(candidate, candidateLc);
        }
        if constexpr (isMcRec) {
          // isSignal =
          //   (TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_sigmac::DecayType::Sc0ToPKPiPi) && chargeCand == 0) ||
          //   (TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_sigmac::DecayType::ScplusplusToPKPiPi) && std::abs(chargeCand) == 2);
          isSignal =
            (std::abs(candidate.flagMcMatchRec()) == (1 << aod::hf_cand_sigmac::DecayType::Sc0ToPKPiPi) && chargeCand == chargeScZero) ||
            (std::abs(candidate.flagMcMatchRec()) == (1 << aod::hf_cand_sigmac::DecayType::ScplusplusToPKPiPi) && std::abs(chargeCand) == chargeScPlusPlus);

          auto trackPos1 = candidateLc.template prong0_as<aod::TracksWMc>();
          auto trackPos2 = candidateLc.template prong2_as<aod::TracksWMc>();
          if (calTrkEff && countCand == 1 && (isSignal || !calEffEventWithCand)) {
            calculateTrkEff(trackPos1, trackPos2, *mcParticles);
          }
          registry.fill(HIST("hPtProng1"), candidate.template prong1_as<aod::TracksWMc>().pt());
        } else {
          registry.fill(HIST("hPtProng1"), candidate.prong1().pt());
        }
        registry.fill(HIST("hPtProng0"), ptCandLc);
      } else {
        selLcPKPi = candidate.isSelLcToPKPi() >= selectionFlagLc;
        selLcPiKP = candidate.isSelLcToPiKP() >= selectionFlagLc;
        if (selLcPKPi) {
          const auto& probs = candidate.mlProbLcToPKPi();
          fillMlOutput(probs, outputMlPKPi);
          massCandPKPi = hfHelper.invMassLcToPKPi(candidate);
        }
        if (selLcPiKP) {
          const auto& probs = candidate.mlProbLcToPiKP();
          fillMlOutput(probs, outputMlPiKP);
          massCandPiKP = hfHelper.invMassLcToPiKP(candidate);
        }
        auto trackPos1 = candidate.template prong0_as<TrackType>();
        auto trackPos2 = candidate.template prong2_as<TrackType>();
        chargeCand = trackPos1.sign();
        if constexpr (isMcRec) {
          isSignal = std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi;
          if (calTrkEff && countCand == 1 && (isSignal || !calEffEventWithCand)) {
            calculateTrkEff(trackPos1, trackPos2, *mcParticles);
          }
        }
        registry.fill(HIST("hPtProng0"), candidate.ptProng0());
        registry.fill(HIST("hPtProng1"), candidate.ptProng1());
        registry.fill(HIST("hPtProng2"), candidate.ptProng2());
      }

      if (isSignal) {
        registry.fill(HIST("hPtCandSig"), ptCand);
        registry.fill(HIST("hEtaSig"), ptCand);
        registry.fill(HIST("hPhiSig"), phiCand);
        registry.fill(HIST("hYSig"), yCand);
      }

      if (selLcPKPi) {
        registry.fill(HIST("hMassLcVsPt"), massCandPKPi, ptCand, efficiencyWeightCand);
        registry.fill(HIST("hMassScVsPtVsSign"), massCandPKPi, ptCand, chargeCand, efficiencyWeightCand);
        registry.fill(HIST("hMassLcData"), massCandPKPi, efficiencyWeightCand);
        registry.fill(HIST("hSelectionStatusLcToPKPi"), selLcPKPi);
        if (isPrompt) {
          registry.fill(HIST("hPtCandSigPrompt"), ptCand);
          registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), ptCand, multiplicityFT0M);
        } else if (isNonPrompt) {
          registry.fill(HIST("hPtCandSigNonPrompt"), ptCand);
          registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), ptCand, multiplicityFT0M);
        }

        entryCandCandRecoInfo(massCandPKPi, ptCand, outputMlPKPi[0], outputMlPKPi[1]);
        entryCandCandGenInfo(isPrompt);
        if (!skipMixedEventTableFilling) {
          entryCand(candidate.phi(), etaCand, ptCand, massCandPKPi, poolBin, gCollisionId, timeStamp);
          entryCandCharge(chargeCand);
        }
      }

      if (selLcPiKP) {
        registry.fill(HIST("hMassLcVsPt"), massCandPiKP, ptCand, efficiencyWeightCand);
        registry.fill(HIST("hMassScVsPtVsSign"), massCandPKPi, ptCand, chargeCand, efficiencyWeightCand);
        registry.fill(HIST("hMassLcData"), massCandPiKP, efficiencyWeightCand);
        registry.fill(HIST("hSelectionStatusLcToPiKP"), selLcPiKP);
        if (isPrompt) {
          registry.fill(HIST("hPtCandSigPrompt"), ptCand);
          registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), ptCand, multiplicityFT0M);
        } else if (isNonPrompt) {
          registry.fill(HIST("hPtCandSigNonPrompt"), ptCand);
          registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), ptCand, multiplicityFT0M);
        }
        entryCandCandRecoInfo(massCandPiKP, ptCand, outputMlPiKP[0], outputMlPiKP[1]);
        entryCandCandGenInfo(isPrompt);
        if (!skipMixedEventTableFilling) {
          entryCand(candidate.phi(), etaCand, ptCand, massCandPiKP, poolBin, gCollisionId, timeStamp);
          entryCandCharge(chargeCand);
        }
      }

      registry.fill(HIST("hCandBin"), poolBin);
      // Correlation with hadrons
      for (const auto& track : tracks) {
        // Remove Lc daughters by checking track indices
        if constexpr (!isCandSc) {
          if ((candidate.prong0Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex()) || (candidate.prong2Id() == track.globalIndex())) {
            if (!storeAutoCorrelationFlag) {
              continue;
            }
            correlationStatus = true;
          }
        } else {
          const auto& candidateLc = candidate.template prongLc_as<CandsLcData>();
          if ((candidateLc.prong0Id() == track.globalIndex()) || (candidateLc.prong1Id() == track.globalIndex()) || (candidateLc.prong2Id() == track.globalIndex()) || (candidate.prong1Id() == track.globalIndex())) {
            if (!storeAutoCorrelationFlag) {
              continue;
            }
            correlationStatus = true;
          }
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
        if constexpr (isMcRec) {
          if (calTrkEff && countCand == 1 && (isSignal || !calEffEventWithCand) && track.has_mcParticle()) {
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            if (!mcParticle.isPhysicalPrimary() && isRecTrkPhyPrimary)
              continue;

            auto motherTrk = mcParticles->iteratorAt(mcParticle.mothersIds()[0]);
            if (std::abs(motherTrk.pdgCode()) == kLambdaCPlus)
              continue;

            registry.fill(HIST("hPtTracksVsSignRec"), track.pt(), track.sign());
            if (std::abs(mcParticle.pdgCode()) == kProton)
              registry.fill(HIST("hPtTracksVsSignRecTrue"), track.pt(), track.sign());
          }
        }

        if (selLcPKPi) {
          fillCorrelationTable<isMcRec>(fillTrkPID, track, candidate, outputMlPKPi, poolBin, correlationStatus, yCand, chargeCand, *mcParticles);
          entryCandHadronRecoInfo(massCandPKPi, false);
        }
        if (selLcPiKP) {
          fillCorrelationTable<isMcRec>(fillTrkPID, track, candidate, outputMlPiKP, poolBin, correlationStatus, yCand, chargeCand, *mcParticles);
          entryCandHadronRecoInfo(massCandPiKP, false);
        }

        if (countCand == 1) {
          if (!skipMixedEventTableFilling) {
            entryHadron(track.phi(), track.eta(), track.pt() * track.sign(), poolBin, gCollisionId, timeStamp);
            if (fillTrkPID) {
              entryTrkPID(track.tpcNSigmaPr(), track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tofNSigmaPr(), track.tofNSigmaKa(), track.tofNSigmaPi());
            }
            registry.fill(HIST("hTracksBin"), poolBin);
          }
        }
      } // end Hadron Tracks loop
      countCand++;
    } // end outer Lc loop
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), multiplicityFT0M);
  }

  template <bool isMcRec, bool isCandSc, typename CollisionType, typename CandType, typename TrackType>
  void doMixEvent(CollisionType const& collisions,
                  TrackType const& tracks,
                  CandType const& candidates,
                  aod::McParticles const* mcParticles = nullptr)
  {

    if (candidates.size() == 0) {
      return;
    }

    double yCand = -999.;
    double ptCand = -999.;
    int8_t chargeCand = 3;
    double massCandPKPi = -999.0;
    double massCandPiKP = -999.0;
    bool selLcPKPi = false;
    bool selLcPiKP = false;
    const int chargeScPlusPlus = 2;
    const int chargeScZero = 0;

    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<CollisionType, CandType, TrackType, BinningType> pairData{corrBinning, numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      poolBinLc = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M()));
      registry.fill(HIST("hMultFT0M"), c1.multFT0M());
      registry.fill(HIST("hZvtx"), c1.posZ());
      registry.fill(HIST("hTracksPoolBin"), poolBin);
      registry.fill(HIST("hLcPoolBin"), poolBinLc);
      for (const auto& [candidate, assocParticle] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        yCand = estimateY<isCandSc>(candidate);
        ptCand = candidate.pt();
        if constexpr (isMcRec) {
          isPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
          isNonPrompt = candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
        }

        if constexpr (isCandSc) {
          const auto& candidateLc = candidate.template prongLc_as<CandsLcData>();
          chargeCand = candidate.charge();
          selLcPKPi = (candidateLc.isSelLcToPKPi() >= selectionFlagLc) && (candidate.statusSpreadLcMinvPKPiFromPDG());
          selLcPiKP = (candidateLc.isSelLcToPiKP() >= selectionFlagLc) && (candidate.statusSpreadLcMinvPiKPFromPDG());
          if (selLcPKPi) {
            const auto& probs = candidateLc.mlProbLcToPKPi();
            fillMlOutput(probs, outputMlPKPi);
            massCandPKPi = hfHelper.invMassScRecoLcToPKPi(candidate, candidateLc);
          }
          if (selLcPiKP) {
            const auto& probs = candidateLc.mlProbLcToPiKP();
            fillMlOutput(probs, outputMlPiKP);
            massCandPiKP = hfHelper.invMassScRecoLcToPKPi(candidate, candidateLc);
          }
          if constexpr (isMcRec) {
            isSignal =
              (TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_sigmac::DecayType::Sc0ToPKPiPi) && chargeCand == chargeScZero) ||
              (TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_sigmac::DecayType::ScplusplusToPKPiPi) && std::abs(chargeCand) == chargeScPlusPlus);
          }
        } else {
          selLcPKPi = candidate.isSelLcToPKPi() >= selectionFlagLc;
          selLcPiKP = candidate.isSelLcToPiKP() >= selectionFlagLc;
          if (selLcPKPi) {
            const auto& probs = candidate.mlProbLcToPKPi();
            fillMlOutput(probs, outputMlPKPi);
            massCandPKPi = hfHelper.invMassLcToPKPi(candidate);
          }
          if (selLcPiKP) {
            const auto& probs = candidate.mlProbLcToPiKP();
            fillMlOutput(probs, outputMlPiKP);
            massCandPiKP = hfHelper.invMassLcToPiKP(candidate);
          }
          auto trackPos1 = candidate.template prong0_as<TrackType>();
          chargeCand = trackPos1.sign();
          if constexpr (isMcRec) {
            isSignal = std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi;
          }
        }

        if (!assocParticle.isGlobalTrackWoDCA() || std::abs(yCand) > yCandMax) {
          continue;
        }

        if (pidTrkApplied) {
          if (!passPIDSelection(assocParticle, trkPIDspecies, pidTPCMax, pidTOFMax, tofPIDThreshold, forceTOF))
            continue;
        }

        if (selLcPKPi) {
          fillCorrelationTable<isMcRec>(fillTrkPID, assocParticle, candidate, outputMlPKPi, poolBin, correlationStatus, yCand, chargeCand, *mcParticles);
          entryCandHadronRecoInfo(massCandPKPi, false);

          if (isPrompt) {
            registry.fill(HIST("hPtCandMcRecSigPrompt"), ptCand);
            registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), ptCand, 0);
          } else if (isNonPrompt) {
            registry.fill(HIST("hPtCandMcRecSigNonPrompt"), ptCand);
            registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), ptCand, 0);
          }
        }

        if (selLcPiKP) {
          fillCorrelationTable<isMcRec>(fillTrkPID, assocParticle, candidate, outputMlPiKP, poolBin, correlationStatus, yCand, chargeCand, *mcParticles);
          entryCandHadronRecoInfo(massCandPiKP, false);

          if (isPrompt) {
            registry.fill(HIST("hPtCandMcRecSigPrompt"), ptCand);
            registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), ptCand, 0);
          } else if (isNonPrompt) {
            registry.fill(HIST("hPtCandMcRecSigNonPrompt"), ptCand);
            registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), ptCand, 0);
          }
        }
      }
    }
  }

  template <bool isCandSc, typename CollisionType, typename PartType>
  void doSameEventMcGen(CollisionType const& mcCollision, PartType const& mcParticles)
  {

    int counterCharmCand = 0;
    registry.fill(HIST("hMcEvtCount"), 0);
    BinningTypeMcGen corrBinningMcGen{{binsZVtx, binsMultiplicityMc}, true};
    poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), mcCollision.multMCFT0A()));
    registry.fill(HIST("hMultFT0AMcGen"), mcCollision.multMCFT0A());

    // find leading particle
    if (correlateLcWithLeadingParticle) {
      leadingIndex = findLeadingParticleMcGen(mcParticles, etaTrackMax.value, ptTrackMin.value);
    }
    // Mc Gen level
    for (const auto& particle : mcParticles) {

      double massCand = -999.0;
      bool isCandFound = isCandSc ? matchCandAndMass<true>(particle, massCand) : matchCandAndMass<false>(particle, massCand);
      if (!isCandFound) {
        continue;
      }
      double yCand = RecoDecay::y(particle.pVector(), massCand);

      if (std::abs(yCand) > yCandGenMax || particle.pt() < ptCandMin) {
        continue;
      }
      registry.fill(HIST("hCandBin"), poolBin);
      registry.fill(HIST("hPtCandMcGen"), particle.pt());
      registry.fill(HIST("hEtaMcGen"), particle.eta());
      registry.fill(HIST("hPhiMcGen"), RecoDecay::constrainAngle(particle.phi(), -PIHalf));
      registry.fill(HIST("hYMcGen"), yCand);

      int8_t chargeLc = pdg->GetParticle(particle.pdgCode())->Charge(); // Retrieve charge
      if (chargeLc != 0) {
        chargeLc = chargeLc / std::abs(chargeLc);
      }

      isPrompt = particle.originMcGen() == RecoDecay::OriginType::Prompt;
      isNonPrompt = particle.originMcGen() == RecoDecay::OriginType::NonPrompt;
      if (isPrompt) {
        registry.fill(HIST("hPtCandMcGenPrompt"), particle.pt());
        registry.fill(HIST("hPtCandVsChargeMcGenPrompt"), particle.pt(), chargeLc);
      } else if (isNonPrompt) {
        registry.fill(HIST("hPtCandMcGenNonPrompt"), particle.pt());
        registry.fill(HIST("hPtCandVsChargeMcGenNonPrompt"), particle.pt(), chargeLc);
      }

      static constexpr std::size_t NDaughtersSc{4u};
      static constexpr std::size_t NDaughtersLc{3u};
      std::vector<int> listDaughters{};
      listDaughters.clear();
      const std::size_t nDaughtersExpected = isCandSc ? NDaughtersSc : NDaughtersLc;

      if (isCandSc) {
        if (massCand == o2::constants::physics::MassSigmaC0 || massCand == o2::constants::physics::MassSigmaCStar0) {
          std::array<int, NDaughtersSc> arrDaughSc0PDG = {kProton, -kKPlus, kPiPlus, kPiMinus};
          RecoDecay::getDaughters(particle, &listDaughters, arrDaughSc0PDG, 2);
        } else {
          std::array<int, NDaughtersSc> arrDaughScPlusPDG = {kProton, -kKPlus, kPiPlus, kPiPlus};
          RecoDecay::getDaughters(particle, &listDaughters, arrDaughScPlusPDG, 2);
        }
      } else {
        std::array<int, NDaughtersLc> arrDaughLcPDG = {kProton, -kKPlus, kPiPlus};
        RecoDecay::getDaughters(particle, &listDaughters, arrDaughLcPDG, 2);
      }

      int counterDaughters = 0;
      std::vector<int> prongsId(nDaughtersExpected);
      if (listDaughters.size() == nDaughtersExpected) {
        for (const auto& dauIdx : listDaughters) {
          auto daughI = mcParticles.rawIteratorAt(dauIdx - mcParticles.offset());
          counterDaughters += 1;
          prongsId[counterDaughters - 1] = daughI.globalIndex();
        }
      }
      counterCharmCand++;

      // Lc Hadron correlation dedicated section
      // if it's a Lc particle, search for Hadron and evalutate correlations
      registry.fill(HIST("hcountCandtriggersMcGen"), 0, particle.pt()); // to count trigger Lc for normalisation
      for (const auto& particleAssoc : mcParticles) {
        if (std::abs(particleAssoc.eta()) > etaTrackMax || particleAssoc.pt() < ptTrackMin || particleAssoc.pt() > ptTrackMax) {
          continue;
        }

        if (std::find(prongsId.begin(), prongsId.end(), particleAssoc.globalIndex()) != prongsId.end()) {
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

        int trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, particleAssoc, true);
        int8_t chargeAssoc = pdg->GetParticle(particleAssoc.pdgCode())->Charge(); // Retrieve charge
        chargeAssoc = chargeAssoc / std::abs(chargeAssoc);
        registry.fill(HIST("hPtParticleAssocMcGen"), particleAssoc.pt());
        entryCandHadronPair(getDeltaPhi(particleAssoc.phi(), particle.phi()),
                            particleAssoc.eta() - particle.eta(),
                            particle.pt(),
                            particleAssoc.pt() * chargeAssoc,
                            poolBin,
                            correlationStatus);
        entryCandHadronPairY(particleAssoc.y() - yCand);
        entryCandHadronRecoInfo(massCand, true);
        entryCandHadronGenInfo(isPrompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
        entryPairCandCharge(chargeLc);
      } // end inner loop
    } // end outer loop
    registry.fill(HIST("hcountCandHadronPerEvent"), counterCharmCand);
    registry.fill(HIST("hZvtx"), mcCollision.posZ());
  }

  //}

  /// Lc-hadron correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processDataLc(SelCollisions::iterator const& collision,
                     TracksData const& tracks,
                     CandsLcDataFiltered const& candidates,
                     aod::BCsWithTimestamps const&)
  {
    doSameEvent<false, false>(collision, tracks, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processDataLc, "Process data", true);

  void processDataSc(SelCollisions::iterator const& collision,
                     TracksData const& tracks,
                     aod::Tracks const&,
                     aod::HfCandSc const& candidates,
                     CandsLcData const&,
                     aod::BCsWithTimestamps const&) // MUST be last among index-compatible
  {
    doSameEvent<false, true>(collision, tracks, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processDataSc, "Process data Sc", false);

  /// Lc-Hadron correlation process starts for McRec
  void processMcRecLc(SelCollisions::iterator const& collision,
                      TracksWithMc const& tracks,
                      CandsLcMcRecFiltered const& candidates,
                      aod::McParticles const& mcParticles)
  {
    doSameEvent<true, false>(collision, tracks, candidates, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processMcRecLc, "Process Mc Reco mode", false);

  /// Lc-Hadron correlation process starts for McRec
  void processMcRecSc(SelCollisions::iterator const& collision,
                      TracksWithMc const& tracks,
                      aod::TracksWMc const&,
                      CandsScMcRec const& candidates,
                      CandsLcData const&,
                      aod::McParticles const& mcParticles)
  {
    doSameEvent<true, true>(collision, tracks, candidates, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processMcRecSc, "Process Mc Reco mode", false);

  void processDataMixedEventSc(SelCollisions const& collisions,
                               TracksData const& tracks,
                               aod::Tracks const&,
                               aod::HfCandSc const& candidates,
                               CandsLcData const&)
  {
    doMixEvent<false, true>(collisions, tracks, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processDataMixedEventSc, "Process Mixed Event Data", false);

  void processDataMixedEventLc(SelCollisions const& collisions,
                               CandsLcDataFiltered const& candidates,
                               TracksData const& tracks)
  {

    doMixEvent<false, false>(collisions, tracks, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processDataMixedEventLc, "Process Mixed Event Data", false);

  void processMcRecMixedEventSc(SelCollisions const& collisions,
                                TracksWithMc const& tracks,
                                aod::TracksWMc const&,
                                soa::Join<aod::HfCandSc,
                                          aod::HfCandScMcRec> const& candidates,
                                CandsLcData const&,
                                aod::McParticles const& mcParticles)
  {
    doMixEvent<true, true>(collisions, tracks, candidates, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processMcRecMixedEventSc, "Process Mixed Event McRec", false);

  void processMcRecMixedEventLc(SelCollisions const& collisions,
                                CandsLcMcRecFiltered const& candidates,
                                TracksWithMc const& tracks,
                                aod::McParticles const& mcParticles)
  {
    doMixEvent<true, false>(collisions, tracks, candidates, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processMcRecMixedEventLc, "Process Mixed Event McRec", false)

  /// Lc-Hadron correlation pair builder - for Mc Gen-level analysis
  void processMcGenLc(SelCollisionsMc::iterator const& mcCollision,
                      CandidatesLcMcGen const& mcParticles)
  {
    doSameEventMcGen<false>(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processMcGenLc, "Process Mc Gen Lc mode", false);

  void processMcGenSc(SelCollisionsMc::iterator const& mcCollision,
                      CandidatesScMcGen const& mcParticles)
  {
    doSameEventMcGen<true>(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processMcGenSc, "Process Mc Gen Sc mode", false);

  void processMcGenMixedEvent(SelCollisionsMc const& collisions,
                              CandidatesLcMcGen const& mcParticles)
  {
    BinningTypeMcGen corrBinningMcGen{{binsZVtx, binsMultiplicityMc}, true};
    auto tracksTuple = std::make_tuple(mcParticles, mcParticles);
    Pair<SelCollisionsMc, CandidatesLcMcGen, CandidatesLcMcGen, BinningTypeMcGen> pairMcGen{corrBinningMcGen, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcGen) {
      poolBin = corrBinningMcGen.getBin(std::make_tuple(c1.posZ(), c1.multMCFT0A()));
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
        bool isPrompt = candidate.originMcGen() == RecoDecay::OriginType::Prompt;
        entryCandHadronPair(getDeltaPhi(particleAssoc.phi(), candidate.phi()),
                            particleAssoc.eta() - candidate.eta(),
                            candidate.pt() * chargeLc / std::abs(chargeLc),
                            particleAssoc.pt() * chargeAssoc / std::abs(chargeAssoc),
                            poolBin,
                            correlationStatus);
        entryCandHadronPairY(particleAssoc.y() - yL);
        entryCandHadronRecoInfo(MassLambdaCPlus, true);
        entryCandHadronGenInfo(isPrompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorLcScHadrons, processMcGenMixedEvent, "Process Mixed Event McGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorLcScHadronsSelection>(cfgc),
                      adaptAnalysisTask<HfCorrelatorLcScHadrons>(cfgc)};
}
