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

/// \file correlatorXicHadrons.cxx
/// \brief Xic-Hadrons correlator task - data-like, Mc-Reco and Mc-Gen analyses
/// \author Ravindra Singh <ravindra.singh@cern.ch>

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/DecayChannelsLegacy.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFC/Utils/utilsCorrelations.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

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

// ============================================================================
// ENUMS FOR MAGIC NUMBERS - REPLACES MAGIC NUMBERS WITH TYPED CONSTANTS
// ============================================================================

/// Enum for candidate types (Xic0 vs XicPlus)
enum class CandidateTypeXic : int8_t {
  Xic0 = 0,
  XicPlus = 1
};

/// Enum for particle types (particle vs antiparticle)
enum class ParticleType : int8_t {
  Particle = 1,
  AntiParticle = -1
};

/// Enum for V0 lambda types
enum class V0LambdaType : int8_t {
  Lambda = 1,
  AntiLambda = -1
};

/// Enum for decay daughters count
enum class XicDecayDaughtersCount : size_t {
  Xic0DaughtersCount = 4u,
  XicPlusDaughtersCount = 5u
};

/// Enum for PDG charge scale factor
enum class PDGChargeScale : size_t {
  Scale = 3u
};

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

///
/// Returns deltaPhi values in range [-pi/2., 3.*pi/2.]
///
double getDeltaPhi(double phiXic, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiXic, -PIHalf);
}

// definition of ME variables
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
using BinningTypeMcGen = ColumnBinningPolicy<aod::mccollision::PosZ, o2::aod::mult::MultMCFT0A>;

// ============================================================================
// SELECTION STRUCT - FOR COLLISION SELECTION WITH XIC
// ============================================================================

struct HfCorrelatorXicHadronsSelection {
  Produces<aod::LcSelection> candSel; // using LcSelection table to avoid duplication for similair table, name will be changed later

  Configurable<bool> useSel8{"useSel8", true, "Flag for applying sel8 for collision selection"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<bool> doSelXicCollision{"doSelXicCollision", true, "Select collisions with at least one Xic"};
  Configurable<int> selectionFlagXic{"selectionFlagXic", 1, "Selection Flag for Xic"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};

  struct : ConfigurableGroup {
    Configurable<float> cfgV0radiusMin{"cfgV0radiusMin", 1.2, "minimum decay radius"};
    Configurable<float> cfgDCAPosToPVMin{"cfgDCAPosToPVMin", 0.05, "minimum DCA to PV for positive track"};
    Configurable<float> cfgDCANegToPVMin{"cfgDCANegToPVMin", 0.2, "minimum DCA to PV for negative track"};
    Configurable<float> cfgV0CosPA{"cfgV0CosPA", 0.995, "minimum v0 cosine"};
    Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};
    Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
    Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};
    Configurable<float> cfgPV{"cfgPV", 10., "maximum z-vertex"};
    Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
    Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  } cfgV0;

  SliceCache cache;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  // Xic+ tables
  using CandsXicPlusDataFiltered = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi>>;
  using CandsXicPlusMcRecFiltered = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi, aod::HfCandXicMcRec>>;
  // Xic0 tables
  using CandsXic0DataFiltered = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf>>;
  using CandsXic0McRecFiltered = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec>>;

  // MC Gen
  using CandidatesXicPlusMcGen = soa::Join<aod::McParticles, aod::HfCandXicMcGen>;
  using CandidatesXic0McGen = soa::Join<aod::McParticles, aod::HfXicToXiPiMCGen>;

  // filter on selection of Xic and decay channel Xic->PKPi
  // Filter xicPlusFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << o2::aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_xic::isSelXicToXiPiPi || selectionFlagXic <= 0);
  // Filter xicPlusFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << o2::aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= selectionFlagXic);
  // Filter xic0Filter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_toxipi::resultSelections >= selectionFlagXic);
  Filter xicPlusFilter = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= selectionFlagXic;
  Filter xic0Filter = aod::hf_sel_toxipi::resultSelections == true;
  template <bool isXicPlus, typename CollType, typename CandType>
  void selectionCollision(CollType const& collision, CandType const& candidates)
  {
    bool isSelColl = true;
    bool isCandFound = false;
    bool isSel8 = true;
    bool isNosameBunchPileUp = true;
    double yCand = -999.;
    double massCand = -999.;
    double ptCand = -999;
    if (doSelXicCollision) {
      for (const auto& candidate : candidates) {
        // For both XicPlus and Xic0
        if constexpr (isXicPlus) {
          massCand = o2::constants::physics::MassXiCPlus;
          ptCand = candidate.pt();
          yCand = candidate.y(massCand);
        } else {
          massCand = o2::constants::physics::MassXiC0;
          ptCand = RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
          yCand = candidate.kfRapXic();
        }

        if (std::abs(yCand) > yCandMax || ptCand < ptCandMin) {
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

  template <typename CandType>
  void selectionCollisionMcGen(CandType const& mcParticles)
  {
    bool isCandFound = false;
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != kXiCPlus && std::abs(particle.pdgCode()) != kXiC0) {
        continue;
      }

      double const massCand = (std::abs(particle.pdgCode()) == kXiC0) ? o2::constants::physics::MassXiC0 : o2::constants::physics::MassXiCPlus;
      double const yCand = RecoDecay::y(particle.pVector(), massCand);
      if (std::abs(yCand) > yCandMax || particle.pt() < ptCandMin) {
        isCandFound = false;
        continue;
      }

      isCandFound = true;
      break;
    }
    candSel(isCandFound);
  }

  template <typename TCollision, typename V0>
  bool selectionV0(TCollision const& collision, V0 const& candidate)
  {
    if (candidate.v0radius() < cfgV0.cfgV0radiusMin) {
      return false;
    }
    if (std::abs(candidate.dcapostopv()) < cfgV0.cfgDCAPosToPVMin) {
      return false;
    }
    if (std::abs(candidate.dcanegtopv()) < cfgV0.cfgDCANegToPVMin) {
      return false;
    }
    if (candidate.v0cosPA() < cfgV0.cfgV0CosPA) {
      return false;
    }
    if (std::abs(candidate.dcaV0daughters()) > cfgV0.cfgDCAV0Dau) {
      return false;
    }
    if (candidate.pt() < cfgV0.cfgV0PtMin) {
      return false;
    }
    if (std::abs(candidate.yLambda()) > yCandMax) {
      return false;
    }
    if (candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda > cfgV0.cfgV0LifeTime) {
      return false;
    }

    return true;
  }

  template <typename TCollision>
  bool eventSelV0(TCollision collision)
  {
    if (!collision.sel8()) {
      return 0;
    }

    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (std::abs(collision.posZ()) > cfgV0.cfgPV) {
      return 0;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    if (collision.trackOccupancyInTimeRange() > cfgV0.cfgMaxOccupancy || collision.trackOccupancyInTimeRange() < cfgV0.cfgMinOccupancy) {
      return 0;
    }

    return 1;
  } // event selection V0

  void processV0Selection(SelCollisions::iterator const& collision,
                          aod::V0Datas const& V0s)
  {
    bool isCandFound = false;
    const int64_t kMinV0Candidates = 1;

    if (!eventSelV0(collision)) {
      candSel(isCandFound);
      return;
    }
    if (V0s.size() < kMinV0Candidates) {
      candSel(isCandFound);
      return;
    }
    for (const auto& v0 : V0s) {
      if (selectionV0(collision, v0)) {
        isCandFound = true;
        break;
      }
    }
    candSel(isCandFound);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadronsSelection, processV0Selection, "Process V0 Collision Selection for Data", false);

  void processXicPlusSelection(SelCollisions::iterator const& collision,
                               CandsXicPlusDataFiltered const& candidates)
  {
    selectionCollision<true>(collision, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadronsSelection, processXicPlusSelection, "Process XicPlus Collision Selection for Data", true);

  void processXic0Selection(SelCollisions::iterator const& collision,
                            CandsXic0DataFiltered const& candidates)
  {
    selectionCollision<false>(collision, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadronsSelection, processXic0Selection, "Process Xic0 Collision Selection for Data", false);

  void processXicPlusSelectionMcRec(SelCollisions::iterator const& collision,
                                    CandsXicPlusMcRecFiltered const& candidates)
  {
    selectionCollision<true>(collision, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadronsSelection, processXicPlusSelectionMcRec, "Process XicPlus Selection McRec", false);

  void processXic0SelectionMcRec(SelCollisions::iterator const& collision,
                                 CandsXic0McRecFiltered const& candidates)
  {
    selectionCollision<false>(collision, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadronsSelection, processXic0SelectionMcRec, "Process Xic0 Selection McRec", false);

  void processXicPlusSelectionMcGen(aod::McCollision const&,
                                    CandidatesXicPlusMcGen const& mcParticles)
  {
    selectionCollisionMcGen(mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadronsSelection, processXicPlusSelectionMcGen, "Process XicPlus Selection McGen", false);

  void processXic0SelectionMcGen(aod::McCollision const&,
                                 CandidatesXic0McGen const& mcParticles)
  {
    selectionCollisionMcGen(mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadronsSelection, processXic0SelectionMcGen, "Process Xic0 Selection McGen", false);
};

// Xic-Hadron correlation pair builder
struct HfCorrelatorXicHadrons {
  // using Lc correlation table to avoid duplication for similair tables, name will be changed later
  Produces<aod::LcHadronPair> entryCandHadronPair;
  Produces<aod::LcHadronPairY> entryCandHadronPairY;
  Produces<aod::LcHadronPairTrkPID> entryCandHadronPairTrkPID;
  Produces<aod::LcHadronRecoInfo> entryCandHadronRecoInfo;
  Produces<aod::LcHadronMlInfo> entryCandHadronMlInfo;
  Produces<aod::LcHadronGenInfo> entryCandHadronGenInfo;
  Produces<aod::LcRecoInfo> entryCandRecoInfo;
  Produces<aod::LcGenInfo> entryCandCandGenInfo;
  Produces<aod::TrkRecInfoLc> entryTrackRecoInfo;
  Produces<aod::Lc> entryCand;
  Produces<aod::Hadron> entryHadron;
  Produces<aod::LcHadronTrkPID> entryTrkPID;
  Produces<aod::CandChargePair> entryPairCandCharge;
  Produces<aod::CandCharge> entryCandCharge;
  Produces<aod::CandHadronInvMass> entryXicHadronInvMass;
  Produces<aod::PairedV0InvMass> entryPairedV0InvMass;
  Produces<aod::V0InvMass> entryV0InvMass;

  struct : ConfigurableGroup {
    Configurable<int> selectionFlagXic{"selectionFlagXic", 1, "Selection flag for Xic"};
    Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "number of events mixed in ME process"};
    Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying Xic efficiency weights"};
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
    Configurable<std::vector<double>> binsPtXic{"binsPtXic", std::vector<double>{0.0, 2.0, 4.0, 6.0, 8.0, 12.0, 50.0}, "pT bin limits for XicPlus mass plots"};
    Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle"};
    Configurable<std::vector<double>> binsPtEfficiencyXic{"binsPtEfficiencyXic", std::vector<double>{0.0, 2.0, 4.0, 6.0, 8.0, 12.0, 50.0}, "pT bin limits for efficiency"};
    Configurable<std::vector<double>> efficiencyXic{"efficiencyXic", {1., 1., 1., 1., 1., 1.}, "efficiency values for Xic"};
    Configurable<bool> storeAutoCorrelationFlag{"storeAutoCorrelationFlag", false, "Store flag that indicates if the track is paired to its Xic mother instead of skipping it"};
    Configurable<bool> correlateXicWithLeadingParticle{"correlateXicWithLeadingParticle", false, "Switch for correlation of Xic baryons with leading particle only"};
    Configurable<bool> pidTrkApplied{"pidTrkApplied", false, "Apply PID selection for associated tracks"};
    Configurable<std::vector<int>> trkPIDspecies{"trkPIDspecies", std::vector<int>{o2::track::PID::Proton, o2::track::PID::Pion, o2::track::PID::Kaon}, "Trk sel: Particles species for PID"};
    Configurable<std::vector<float>> pidTPCMax{"pidTPCMax", std::vector<float>{3., 0., 0.}, "maximum nSigma TPC"};
    Configurable<std::vector<float>> pidTOFMax{"pidTOFMax", std::vector<float>{3., 0., 0.}, "maximum nSigma TOF"};
    Configurable<float> tofPIDThreshold{"tofPIDThreshold", 0.75, "minimum pT after which TOF PID is applicable"};
    Configurable<bool> fillTrkPID{"fillTrkPID", false, "fill PID information for associated tracks"};
    Configurable<bool> forceTOF{"forceTOF", false, "fill PID information for associated tracks"};
    Configurable<bool> calTrkEff{"calTrkEff", false, "fill histograms to calculate efficiency"};
    Configurable<bool> isRecTrkPhyPrimary{"isRecTrkPhyPrimary", true, "Calculate the efficiency of reconstructed primary physical tracks"};
    Configurable<bool> calEffEventWithCand{"calEffEventWithCand", true, "Calculate the efficiency of Xic candidate"};
    Configurable<float> eventFractionToAnalyze{"eventFractionToAnalyze", -1, "Fraction of events to analyze"};
    Configurable<int> particlePdg{"particlePdg", PDG_t::kProton, "PDG code of particle: kProton(2212), kPiPlus(211), kKPlus(321), kLambda0(3122)"};
    Configurable<int8_t> corrParticle{"corrParticle", 1, "put '0' for physical primary, '1' for indivisual identified particle, '2' for V0s"};
  } cfgXicCand;

  struct : ConfigurableGroup {
    Configurable<float> cfgDaughPrPtMax{"cfgDaughPrPtMax", 5.f, "max. pT daughter proton"};
    Configurable<float> cfgDaughPrPtMin{"cfgDaughPrPtMin", 0.3f, "min. pT daughter proton"};
    Configurable<float> cfgDaughPiPtMax{"cfgDaughPiPtMax", 10.f, "max. pT daughter pion"};
    Configurable<float> cfgDaughPiPtMin{"cfgDaughPiPtMin", 0.3f, "min. pT daughter pion"};
    Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 3.f, "max. TPC nSigma proton"};
    Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 2.f, "max. TPC nSigma pion"};
    Configurable<float> cfgDaughPIDCutsTOFPi{"cfgDaughPIDCutsTOFPi", 2.f, "max. TOF nSigma pion"};
    Configurable<float> cfgHypMassWindow{"cfgHypMassWindow", 0.1f, "single lambda mass selection"};
    Configurable<bool> cfgIsCorrCollMatchV0{"cfgIsCorrCollMatchV0", true, "check if daughter and mother collision are same"};
    Configurable<bool> cfgCalDataDrivenEffPr{"cfgCalDataDrivenEffPr", false, "calculate data driven efficiency of proton using Lambda"};
  } cfgV0;

  SliceCache cache;
  Service<o2::framework::O2DatabasePDG> pdg{};
  int8_t chargeCand = 3;
  int leadingIndex = 0;
  int poolBin = 0;
  int poolBinXic = 0;
  bool correlationStatus = false;
  bool isPrompt = false;
  bool isNonPrompt = false;
  bool isSignal = false;
  TRandom3* rnd = new TRandom3(0);
  std::vector<float> outputMlXic = {-1., -1., -1.};

  // Event Mixing for the Data Mode
  using SelCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::LcSelection>>;
  using SelCollisionsMc = soa::Filtered<soa::Join<aod::McCollisions, aod::LcSelection, aod::MultsExtraMC>>; // collisionFilter applied

  // XicPlus data
  using CandsXicPlusData = soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi, aod::HfMlXicToXiPiPi>;
  using CandsXicPlusDataFiltered = soa::Filtered<CandsXicPlusData>;
  using CandsXicPlusMcRec = soa::Join<aod::HfCandXic, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi, aod::HfMlXicToXiPiPi>;
  using CandsXicPlusMcRecFiltered = soa::Filtered<CandsXicPlusMcRec>;

  // Xic0 data
  using CandsXic0Data = soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfMlToXiPi>;
  using CandsXic0DataFiltered = soa::Filtered<CandsXic0Data>;
  using CandsXic0McRec = soa::Join<aod::HfCandToXiPiKf, aod::HfXicToXiPiMCRec, aod::HfSelToXiPiKf, aod::HfMlToXiPi>;
  using CandsXic0McRecFiltered = soa::Filtered<CandsXic0McRec>;

  // MC Gen
  using CandidatesXicPlusMcGen = soa::Join<aod::McParticles, aod::HfCandXicMcGen>;
  using CandidatesXic0McGen = soa::Join<aod::McParticles, aod::HfXicToXiPiMCGen>;

  using McCollisionsSel = soa::Filtered<soa::Join<aod::McCollisions, aod::LcSelection>>;
  using McParticlesSel = soa::Filtered<aod::McParticles>;

  // Tracks used in Data and MC
  using TracksData = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using TracksWithMc = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  Filter collisionFilter = aod::hf_selection_lc_collision::lcSel == true;
  //  Filter xicPlusFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << o2::aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= cfgXicCand.selectionFlagXic);
  // Filter xicPlusFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << o2::aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_xic::isSelXicToXiPiPi || cfgXicCand.selectionFlagXic <= 0);
  // Filter xic0Filter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_toxipi::resultSelections >= cfgXicCand.selectionFlagXic);
  Filter trackFilter = (nabs(aod::track::eta) < cfgXicCand.etaTrackMax) && (nabs(aod::track::pt) > cfgXicCand.ptTrackMin) && (nabs(aod::track::dcaXY) < cfgXicCand.dcaXYTrackMax) && (nabs(aod::track::dcaZ) < cfgXicCand.dcaZTrackMax);
  Filter xicPlusFilter = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= cfgXicCand.selectionFlagXic;
  Filter xic0Filter = aod::hf_sel_toxipi::resultSelections == true;

  Preslice<aod::McParticles> perTrueCollision = o2::aod::mcparticle::mcCollisionId;
  Preslice<aod::Tracks> perCollisionID = aod::track::collisionId;
  Preslice<aod::HfCandXic> candXicPerCol = aod::hf_cand::collisionId;
  Preslice<aod::HfCandToXiPi> candXic0PerCol = aod::hf_cand_xic0_omegac0::collisionId;

  // configurable axis definition
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0.0f, 2000.0f, 6000.0f, 100000.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis binsZVtx{"binsZVtx", {VARIABLE_WIDTH, -10.0f, -2.5f, 2.5f, 10.0f}, "Mixing bins - z-vertex"};
  ConfigurableAxis binsMultiplicityMc{"binsMultiplicityMc", {VARIABLE_WIDTH, 0.0f, 20.0f, 50.0f, 500.0f}, "Mixing bins - MC multiplicity"};
  ConfigurableAxis binsBdtScore{"binsBdtScore", {100, 0., 1.}, "Bdt output scores"};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 6000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsCandMassXicPlus{"binsCandMassXicPlus", {400, 2., 2.8}, "inv. mass (Xi pi pi) (GeV/#it{c}^{2})"};
  ConfigurableAxis binsCandMassXic0{"binsCandMassXic0", {200, 2.4, 2.8}, "inv. mass (Xi pi) (GeV/#it{c}^{2})"};
  ConfigurableAxis binsNSigmas{"binsNSigmas", {4000, -500., 500.}, "n#sigma"};

  BinningType corrBinning{{binsZVtx, binsMultiplicity}, true};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // XicPlus mass axes
    AxisSpec axisCandMassXicPlus = {binsCandMassXicPlus, "inv. mass (Xi pi pi) (GeV/#it{c}^{2})"};
    // Xic0 mass axes
    AxisSpec axisCandMassXic0 = {binsCandMassXic0, "inv. mass (Xi pi) (GeV/#it{c}^{2})"};

    AxisSpec const axisEta = {binsEta, "#it{eta}"};
    AxisSpec const axisPhi = {binsPhi, "#it{#varphi}"};
    AxisSpec axisPtXic = {static_cast<std::vector<double>>(cfgXicCand.binsPtXic), "#it{p}_{T} Xic (GeV/#it{c})"};
    AxisSpec axisPtHadron = {static_cast<std::vector<double>>(cfgXicCand.binsPtHadron), "#it{p}_{T} Hadron (GeV/#it{c})"};
    AxisSpec axisPtTrack = {500, 0, 50, "#it{p}_{T} Hadron (GeV/#it{c})"};
    AxisSpec const axisMultiplicity = {binsMultiplicity, "Multiplicity"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec const axisPosZ = {binsZVtx, "PosZ"};
    AxisSpec const axisBdtScore = {binsBdtScore, "Bdt score"};
    AxisSpec const axisPoolBin = {binsPoolBin, "PoolBin"};
    AxisSpec const axisRapidity = {100, -2, 2, "Rapidity"};
    AxisSpec const axisNSigma = {binsNSigmas, "n#sigma"};
    AxisSpec axisSign = {5, -2.5, 2.5, "Sign"};
    AxisSpec axisPtV0 = {500, 0., 50.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisMassV0 = {300, 1.05f, 1.2f, "inv. mass (p #pi) (GeV/#it{c}^{2})"};

    // XicPlus histograms
    registry.add("hMassXicPlusVsPt", "XicPlus candidates;inv. mass (Xi pi pi) (GeV/#it{c}^{2});#it{p}_{T};entries", {HistType::kTH2F, {{axisCandMassXicPlus}, {axisPtXic}}});
    registry.add("hMassXicPlusData", "XicPlus candidates;inv. mass (Xi pi pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{axisCandMassXicPlus}}});

    // Xic0 histograms
    registry.add("hMassXic0VsPt", "Xic0 candidates;inv. mass (Xi pi) (GeV/#it{c}^{2});#it{p}_{T};entries", {HistType::kTH2F, {{axisCandMassXic0}, {axisPtXic}}});
    registry.add("hMassXic0Data", "Xic0 candidates;inv. mass (Xi pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{axisCandMassXic0}}});

    // Common histograms
    registry.add("hPtCandXic", "Xic0 candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtXic}});
    registry.add("hEta", "Xic candidates;candidate #it{#eta};entries", {HistType::kTH1F, {axisEta}});
    registry.add("hPhi", "Xic candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {axisPhi}});
    registry.add("hY", "Xic candidates;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hMultFT0M", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("hZvtx", "z vertex;z vertex;entries", {HistType::kTH1F, {{200, -20., 20.}}});
    registry.add("hCandBin", "candidates in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}});
    registry.add("hTracksBin", "Tracks selected in pool Bin;pool Bin;entries", {HistType::kTH1F, {{9, 0., 9.}}});
    // registry.add("hMassXicVsPt", "Xic candidates;inv. mass (Xi pi pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisCandMass}, {axisPtXic}}});
    // registry.add("hMassXicPlusVsPtVsSign", "Xic candidates;inv. mass (#Xi #pi #pi) (GeV/#it{c}^{2});sign;entries", {HistType::kTH3F, {{axisCandMass}, {axisPtXic}, {axisSign}}});
    // registry.add("hMassXicData", "Xic candidates;inv. mass (Xi pi pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{axisCandMass}}});
    registry.add("hXicPoolBin", "Xic candidates pool bin", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hTracksPoolBin", "Particles associated pool bin", {HistType::kTH1F, {axisPoolBin}});
    // Histograms for MC Reco analysis
    registry.add("hMcEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}});
    // registry.add("hMassXicMcRecBkg", "Xic background candidates - MC reco;inv. mass (Xi pi pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{axisCandMass}, {axisPtXic}}});
    registry.add("hPtCandSig", "Xic,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtXic}});
    registry.add("hPtCandSigPrompt", "Xic,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtXic}});
    registry.add("hPtCandSigNonPrompt", "Xic,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtXic}});
    // registry.add("hPtCandMcRecBkg", "Xic,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPtXic}});
    registry.add("hEtaSig", "Xic,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiSig", "Xic,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
    // registry.add("hY", "Xic,Hadron candidates;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hYSig", "Xic,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hPtCandMcRecSigPrompt", "Xic,Hadron candidates Prompt - MC Reco", {HistType::kTH1F, {axisPtXic}});
    registry.add("hPtCandMcRecSigNonPrompt", "Xic,Hadron candidates Non Prompt - MC Reco", {HistType::kTH1F, {axisPtXic}});
    registry.add("hEtaMcRecBkg", "Xic,Hadron candidates - MC Reco", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcRecBkg", "Xic,Hadron candidates - MC Reco", {HistType::kTH1F, {axisPhi}});
    registry.add("hYMcRecBkg", "Xic,Hadron candidates - MC reco;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hFakeTracksMcRec", "Fake tracks - MC Rec", {HistType::kTH1F, {axisPtHadron}});
    registry.add("hPtParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtXic}}});
    registry.add("hPtTracksVsSignRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtTrack}, {axisSign}}});
    registry.add("hPtTracksVsSignRecTrue", "Associated Particle - MC Rec (True)", {HistType::kTH2F, {{axisPtTrack}, {axisSign}}});
    registry.add("hPtTracksVsSignGen", "Associated Particle - MC Gen", {HistType::kTH2F, {{axisPtTrack}, {axisSign}}});
    registry.add("hPtPrimaryParticleAssocVsCandMcRec", "Associated Particle - MC Rec", {HistType::kTH2F, {{axisPtHadron}, {axisPtXic}}});
    registry.add("hPtVsMultiplicityMcRecPrompt", "Multiplicity FT0M - MC Rec Prompt", {HistType::kTH2F, {{axisPtXic}, {axisMultFT0M}}});
    registry.add("hPtVsMultiplicityMcRecNonPrompt", "Multiplicity FT0M - MC Rec Non Prompt", {HistType::kTH2F, {{axisPtXic}, {axisMultFT0M}}});
    // Histograms for MC Gen analysis
    registry.add("hcountCandtriggersMcGen", "Xic trigger particles - MC gen;;N of trigger Xic", {HistType::kTH2F, {{1, -0.5, 0.5}, {axisPtXic}}});
    registry.add("hcountCandHadronPerEvent", "Xic,Hadron particles - MC gen;Number per event;entries", {HistType::kTH1F, {{21, -0.5, 20.5}}});
    registry.add("hPtCandMcGen", "Xic,Hadron particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPtXic}});
    registry.add("hYMcGen", "Xic,Hadron candidates - MC gen;candidate #it{#y};entries", {HistType::kTH1F, {axisRapidity}});
    registry.add("hPtCandMcGenPrompt", "Xic,Hadron particles - MC Gen Prompt", {HistType::kTH1F, {axisPtXic}});
    registry.add("hPtCandVsChargeMcGenPrompt", "Charm Hadron particles - MC Gen Prompt", {HistType::kTH2F, {{axisPtXic}, {axisSign}}});
    registry.add("hPtCandMcGenNonPrompt", "Charm Hadron particles - MC Gen Non Prompt", {HistType::kTH1F, {axisPtXic}});
    registry.add("hPtCandVsChargeMcGenNonPrompt", "Xic,Hadron particles - MC Gen Non Prompt", {HistType::kTH2F, {{axisPtXic}, {axisSign}}});
    registry.add("hPtParticleAssocMcGen", "Associated Particle - MC Gen", {HistType::kTH1F, {axisPtHadron}});
    registry.add("hEtaMcGen", "Xic,Hadron particles - MC Gen", {HistType::kTH1F, {axisEta}});
    registry.add("hPhiMcGen", "Xic,Hadron particles - MC Gen", {HistType::kTH1F, {axisPhi}});
    registry.add("hMultFT0AMcGen", "Xic,Hadron multiplicity FT0A - MC Gen", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("hTOFnSigmaPr", "hTOFnSigmaPr", {HistType::kTH2F, {{axisPtHadron}, {axisNSigma}}});
    registry.add("hTPCnSigmaPr", "hTPCnSigmaPr", {HistType::kTH2F, {{axisPtHadron}, {axisNSigma}}});
    registry.add("hTOFnSigmaPrPiKRej", "hTOFnSigmaPrPiKRej", {HistType::kTH2F, {{axisPtHadron}, {axisNSigma}}});
    registry.add("hTPCnSigmaPrPiKRej", "hTPCnSigmaPrPiKRej", {HistType::kTH2F, {{axisPtHadron}, {axisNSigma}}});

    // Xic-Hadron invariant mass histograms
    registry.add("hXicHadronInvMassVsPt", "Xic+Hadron invariant mass vs combined pT;m_{Xic+h} (GeV/c^2);p_{T,combined} (GeV/c)", {HistType::kTH2F, {{500, 2.5, 5.0}, {100, 0., 50.}}});
    registry.add("hXicLambdaInvMassVsPt", "Xic+Lambda invariant mass vs combined pT;m_{Xic+Lambda} (GeV/c^2);p_{T,combined} (GeV/c)", {HistType::kTH2F, {{500, 3.0, 5.5}, {100, 0., 50.}}});
    registry.add("hXicHadronInvMass", "Xic+Hadron invariant mass;m_{Xic+h} (GeV/c^2);entries", {HistType::kTH1F, {{500, 2.5, 5.0}}});
    registry.add("hXicLambdaInvMass", "Xic+Lambda invariant mass;m_{Xic+Lambda} (GeV/c^2);entries", {HistType::kTH1F, {{500, 3.0, 5.5}}});
    registry.add("hXicHadronPtCombined", "Combined pT (Xic+Hadron);p_{T,combined} (GeV/c);entries", {HistType::kTH1F, {{100, 0., 50.}}});
    registry.add("hXicLambdaPtCombined", "Combined pT (Xic+Lambda);p_{T,combined} (GeV/c);entries", {HistType::kTH1F, {{100, 0., 50.}}});

    // Lambda V0 histograms
    registry.add("hEventLambdaV0", "Lambda, events", {HistType::kTH1F, {{2, 0, 2}}});
    registry.add("hV0Lambda", "V0 Lambda candidates;inv. mass (p #pi) (GeV/#it{c}^{2});GeV/#it{c};GeV/#it{c}", {HistType::kTH3F, {{axisMassV0}, {axisPtV0}, {axisPtHadron}}});
    registry.add("hV0LambdaRefl", "V0 Lambda reflected candidates;inv. mass (p #pi) (GeV/#it{c}^{2});GeV/#it{c};GeV/#it{c}", {HistType::kTH3F, {{axisMassV0}, {axisPtV0}, {axisPtHadron}}});
    registry.add("hV0LambdaPiKRej", "V0 Lambda candidates with #pi K rejection;inv. mass (p #pi) (GeV/#it{c}^{2});GeV/#it{c};GeV/#it{c}", {HistType::kTH3F, {{axisMassV0}, {axisPtV0}, {axisPtHadron}}});
    registry.add("hV0LambdaReflPiKRej", "V0 Lambda reflected candidates with #pi K rejection;inv. mass (p #pi) (GeV/#it{c}^{2});GeV/#it{c};GeV/#it{c}", {HistType::kTH3F, {{axisMassV0}, {axisPtV0}, {axisPtHadron}}});
    registry.add("hV0LambdaMcRec", "McRec V0 Lambda candidates;inv. mass (p #pi) (GeV/#it{c}^{2});GeV/#it{c};GeV/#it{c}", {HistType::kTH3F, {{axisMassV0}, {axisPtV0}, {axisPtHadron}}});
    registry.add("hV0LambdaReflMcRec", "McRec V0 Lambda reflected candidates;inv. mass (p #pi) (GeV/#it{c}^{2});GeV/#it{c};GeV/#it{c}", {HistType::kTH3F, {{axisMassV0}, {axisPtV0}, {axisPtHadron}}});
    registry.add("hV0LambdaPiKRejMcRec", "McRec V0 Lambda candidates with #pi K rejection;inv. mass (p #pi) (GeV/#it{c}^{2});GeV/#it{c};GeV/#it{c}", {HistType::kTH3F, {{axisMassV0}, {axisPtV0}, {axisPtHadron}}});
    registry.add("hV0LambdaReflPiKRejMcRec", "McRec V0 Lambda reflected candidates with #pi K rejection;inv. mass (p #pi) (GeV/#it{c}^{2});GeV/#it{c};GeV/#it{c}", {HistType::kTH3F, {{axisMassV0}, {axisPtV0}, {axisPtHadron}}});

    corrBinning = {{binsZVtx, binsMultiplicity}, true};
  }

  template <typename MlProbType>
  void fillMlOutput(MlProbType const& mlProb, std::vector<float>& outputMl)
  {
    if (mlProb.size() == 0) {
      return;
    }
    for (unsigned int iclass = 0; iclass < cfgXicCand.classMl->size() && iclass < outputMl.size(); iclass++) {
      unsigned int targetIndex = cfgXicCand.classMl->at(iclass);
      if (targetIndex < mlProb.size()) {
        outputMl[iclass] = mlProb[targetIndex];
      }
    }
  }

  //  template <typename CandType>
  //  double estimateY(CandType const& candidate)
  //  {
  //    return HfHelper::yXic(candidate);
  //  }

  float getMassFromPdg(int pdgCode)
  {
    switch (std::abs(pdgCode)) {
      case PDG_t::kProton:
        return MassProton;
      case PDG_t::kPiPlus:
        return MassPiPlus;
      case PDG_t::kKPlus:
        return MassKPlus;
      case PDG_t::kLambda0:
        return MassLambda;
      default:
        return 0.f;
    }
  }

  //  template <typename TrackType>
  //  float getXicType(const auto& candidate)
  //  {
  //    int chargeCand = candidate.sign();
  //    float XicType = 0.f;
  //    if (chargeCand == 0) {
  //      auto bachelorTrack = candidate.template prong0_as<TrackType>();
  //      XicType = (bachelorTrack.sign() > 0) ? 0.5f : -0.5f; // to seprate Xic0 with its anti-particle
  //    } else {
  //      XicType = (chargeCand > 0) ? 1.5f : -1.5f;
  //    }
  //    return XicType;
  //  }

  float getXicTypeMC(const auto& particle)
  {
    int pdgCode = particle.pdgCode();

    switch (pdgCode) {
      case kXiC0:
        return 0.5f; // Xic0
      case -kXiC0:
        return -0.5f; // Xic0-bar
      case kXiCPlus:
        return 1.5f; // XicPlus
      case -kXiCPlus:
        return -1.5f; // XicMinus
      default:
        return 0.f;
    }
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, int pid)
  {
    if (std::abs(pid) == kProton && std::abs(track.tpcNSigmaPr()) > cfgV0.cfgDaughPIDCutsTPCPr) {
      return false;
    }
    if (std::abs(pid) == kPiPlus && (std::abs(track.tpcNSigmaPi()) > cfgV0.cfgDaughPIDCutsTPCPi || std::abs(track.tofNSigmaPi()) > cfgV0.cfgDaughPIDCutsTOFPi)) {
      return false;
    }
    if (std::abs(track.eta()) > cfgXicCand.etaTrackMax) {
      return false;
    }
    if (std::abs(pid) == kProton && track.pt() > cfgV0.cfgDaughPrPtMax) {
      return false;
    }
    if (std::abs(pid) == kProton && track.pt() < cfgV0.cfgDaughPrPtMin) {
      return false;
    }
    if (std::abs(pid) == kPiPlus && track.pt() > cfgV0.cfgDaughPiPtMax) {
      return false;
    }
    if (std::abs(pid) == kPiPlus && track.pt() < cfgV0.cfgDaughPiPtMin) {
      return false;
    }

    return true;
  }

  template <typename assocType>
  float calculateInvMass(float pT1, float eta1, float phi1, assocType const& particle2, float mass1, float mass2)
  {
    ROOT::Math::PtEtaPhiMVector vec1(pT1, eta1, phi1, mass1);
    ROOT::Math::PtEtaPhiMVector vec2(particle2.pt(), particle2.eta(), particle2.phi(), mass2);
    ROOT::Math::PtEtaPhiMVector combined = vec1 + vec2;
    return combined.mass();
  }

  template <typename assocType>
  float calculateCombinedPt(float pT1, float eta1, float phi1, assocType const& particle2)
  {
    ROOT::Math::PtEtaPhiMVector vec1(pT1, eta1, phi1, 0.);
    ROOT::Math::PtEtaPhiMVector vec2(particle2.pt(), particle2.eta(), particle2.phi(), 0.);
    ROOT::Math::PtEtaPhiMVector combined = vec1 + vec2;
    return combined.pt();
  }

  // ============================================================================
  // SAME EVENT WITH V0 LAMBDA PROCESSING
  // ============================================================================

  template <bool IsMcRec, bool isXicPlus, typename CollisionType, typename V0, typename TrackType, typename CandsXics>
  void doSameEventWithV0(CollisionType const& collision, V0 const& v0s, TrackType const& tracks, CandsXics const& candidates, aod::McParticles const* mcParticles = nullptr)
  {
    // Data-driven efficiency calculation for protons using Lambda
    if (cfgV0.cfgCalDataDrivenEffPr) {
      for (const auto& v0 : v0s) {
        auto posTrackV0 = v0.template posTrack_as<TrackType>();
        auto negTrackV0 = v0.template negTrack_as<TrackType>();
        if (cfgV0.cfgIsCorrCollMatchV0 && ((v0.collisionId() != posTrackV0.collisionId()) || (v0.collisionId() != negTrackV0.collisionId()))) {
          continue;
        }

        // Process Lambda (proton + pion)
        if (std::abs(o2::constants::physics::MassLambda - v0.mLambda()) < cfgV0.cfgHypMassWindow) {
          entryHadron(v0.mLambda(), posTrackV0.eta(), posTrackV0.pt() * posTrackV0.sign(), 0, 0, v0.pt());
          entryTrkPID(posTrackV0.tpcNSigmaPr(), posTrackV0.tpcNSigmaKa(), posTrackV0.tpcNSigmaPi(), posTrackV0.tofNSigmaPr(), posTrackV0.tofNSigmaKa(), posTrackV0.tofNSigmaPi());

          if (isSelectedV0Daughter(posTrackV0, kProton) && isSelectedV0Daughter(negTrackV0, kPiPlus)) {
            registry.fill(HIST("hV0Lambda"), v0.mLambda(), v0.pt(), posTrackV0.pt());
            registry.fill(HIST("hV0LambdaRefl"), v0.mAntiLambda(), v0.pt(), negTrackV0.pt());
            registry.fill(HIST("hTPCnSigmaPr"), posTrackV0.pt(), posTrackV0.tpcNSigmaPr());

            if (posTrackV0.hasTOF()) {
              registry.fill(HIST("hTOFnSigmaPr"), posTrackV0.pt(), posTrackV0.tofNSigmaPr());
            }

            if (passPIDSelection(posTrackV0, cfgXicCand.trkPIDspecies, cfgXicCand.pidTPCMax, cfgXicCand.pidTOFMax, cfgXicCand.tofPIDThreshold, cfgXicCand.forceTOF)) {
              registry.fill(HIST("hV0LambdaPiKRej"), v0.mLambda(), v0.pt(), posTrackV0.pt());
              registry.fill(HIST("hV0LambdaReflPiKRej"), v0.mAntiLambda(), v0.pt(), negTrackV0.pt());
              registry.fill(HIST("hTPCnSigmaPrPiKRej"), posTrackV0.pt(), posTrackV0.tpcNSigmaPr());
              if (posTrackV0.hasTOF()) {
                registry.fill(HIST("hTOFnSigmaPrPiKRej"), posTrackV0.pt(), posTrackV0.tofNSigmaPr());
              }
            }
          }
        }

        if (std::abs(o2::constants::physics::MassLambda - v0.mAntiLambda()) < cfgV0.cfgHypMassWindow) {
          entryHadron(v0.mAntiLambda(), negTrackV0.eta(), negTrackV0.pt() * negTrackV0.sign(), 0, 0, v0.pt());
          entryTrkPID(negTrackV0.tpcNSigmaPr(), negTrackV0.tpcNSigmaKa(), negTrackV0.tpcNSigmaPi(), negTrackV0.tofNSigmaPr(), negTrackV0.tofNSigmaKa(), negTrackV0.tofNSigmaPi());

          if (isSelectedV0Daughter(negTrackV0, kProton) && isSelectedV0Daughter(posTrackV0, kPiPlus)) {
            registry.fill(HIST("hV0Lambda"), v0.mAntiLambda(), v0.pt(), negTrackV0.pt());
            registry.fill(HIST("hV0LambdaRefl"), v0.mLambda(), v0.pt(), posTrackV0.pt());
            registry.fill(HIST("hTPCnSigmaPr"), negTrackV0.pt(), negTrackV0.tpcNSigmaPr());

            if (negTrackV0.hasTOF()) {
              registry.fill(HIST("hTOFnSigmaPr"), negTrackV0.pt(), negTrackV0.tofNSigmaPr());
            }

            if (passPIDSelection(negTrackV0, cfgXicCand.trkPIDspecies, cfgXicCand.pidTPCMax, cfgXicCand.pidTOFMax, cfgXicCand.tofPIDThreshold, cfgXicCand.forceTOF)) {
              registry.fill(HIST("hV0LambdaPiKRej"), v0.mAntiLambda(), v0.pt(), negTrackV0.pt());
              registry.fill(HIST("hV0LambdaReflPiKRej"), v0.mLambda(), v0.pt(), posTrackV0.pt());
              registry.fill(HIST("hTPCnSigmaPrPiKRej"), negTrackV0.pt(), negTrackV0.tpcNSigmaPr());
              if (negTrackV0.hasTOF()) {
                registry.fill(HIST("hTOFnSigmaPrPiKRej"), negTrackV0.pt(), negTrackV0.tofNSigmaPr());
              }
            }
          }
        }

        // MC-Reco specific V0 matching
        if constexpr (IsMcRec) {
          if (!v0.has_mcParticle() || !posTrackV0.has_mcParticle() || !negTrackV0.has_mcParticle()) {
            continue;
          }
          auto v0Mc = v0.mcParticle();
          auto posTrack = posTrackV0.mcParticle();
          auto negTrack = negTrackV0.mcParticle();

          if (std::abs(v0Mc.pdgCode()) == kLambda0) {
            if (std::abs(posTrack.pdgCode()) == kProton) {
              registry.fill(HIST("hV0LambdaMcRec"), v0.mLambda(), v0.pt(), posTrackV0.pt());
              registry.fill(HIST("hV0LambdaReflMcRec"), v0.mAntiLambda(), v0.pt(), negTrackV0.pt());

              if (passPIDSelection(posTrackV0, cfgXicCand.trkPIDspecies, cfgXicCand.pidTPCMax, cfgXicCand.pidTOFMax, cfgXicCand.tofPIDThreshold, cfgXicCand.forceTOF)) {
                registry.fill(HIST("hV0LambdaPiKRejMcRec"), v0.mLambda(), v0.pt(), posTrackV0.pt());
                registry.fill(HIST("hV0LambdaReflPiKRejMcRec"), v0.mAntiLambda(), v0.pt(), negTrackV0.pt());
              }
            }
            if (std::abs(negTrack.pdgCode()) == kProton) {
              registry.fill(HIST("hV0LambdaMcRec"), v0.mAntiLambda(), v0.pt(), negTrackV0.pt());
              registry.fill(HIST("hV0LambdaReflMcRec"), v0.mLambda(), v0.pt(), posTrackV0.pt());

              if (passPIDSelection(negTrackV0, cfgXicCand.trkPIDspecies, cfgXicCand.pidTPCMax, cfgXicCand.pidTOFMax, cfgXicCand.tofPIDThreshold, cfgXicCand.forceTOF)) {
                registry.fill(HIST("hV0LambdaPiKRejMcRec"), v0.mAntiLambda(), v0.pt(), negTrackV0.pt());
                registry.fill(HIST("hV0LambdaReflPiKRejMcRec"), v0.mLambda(), v0.pt(), posTrackV0.pt());
              }
            }
          }
        }
      }
    }

    // Main Xic-Lambda correlation loop
    int nTracks = 0;
    float efficiencyWeightCand = 1.;
    int64_t timeStamp = 0;
    bool skipMixedEventTableFilling = false;
    float const multiplicityFT0M = collision.multFT0M();
    int gCollisionId = collision.globalIndex();

    if (candidates.size() == 0) {
      return;
    }

    if (cfgXicCand.eventFractionToAnalyze > 0) {
      if (rnd->Uniform(0, 1) > cfgXicCand.eventFractionToAnalyze) {
        skipMixedEventTableFilling = true;
      }
    }

    if constexpr (!IsMcRec) {
      timeStamp = collision.template bc_as<aod::BCsWithTimestamps>().timestamp();
    }

    poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), multiplicityFT0M));

    // Count good tracks
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (std::abs(track.eta()) > cfgXicCand.etaTrackMax || std::abs(track.dcaXY()) > cfgXicCand.dcaXYTrackMax || std::abs(track.dcaZ()) > cfgXicCand.dcaZTrackMax) {
          continue;
        }
        nTracks++;
      }
    }

    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < cfgXicCand.multMin || nTracks > cfgXicCand.multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    int countCand = 1;

    // Xic-Lambda correlation (same event)
    for (const auto& candidate : candidates) {
      double massCand = -999.0;
      double yCand = -999.0;
      double ptCand = -999.0;
      double absPtCand = -999.0;
      double etaCand = -999.0;
      double phiCand = -999.0;
      bool selXicCand = false;

      // float xicType = getXicType(candidate);

      // Determine mass and rapidity based on Xic type
      if constexpr (!isXicPlus) {
        massCand = candidate.invMassCharmBaryon();
        yCand = candidate.kfRapXic();
        ptCand = -RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon()) * candidate.signDecay();
        etaCand = candidate.etaCharmBaryon();
        phiCand = RecoDecay::phi(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
        selXicCand = candidate.resultSelections() >= cfgXicCand.selectionFlagXic;
        const auto& probs = candidate.mlProbToXiPi();
        fillMlOutput(probs, outputMlXic);
      } else { // XicPlus
        massCand = candidate.invMassXicPlus();
        yCand = candidate.y(o2::constants::physics::MassXiCPlus);
        ptCand = candidate.pt() * candidate.sign();
        etaCand = candidate.eta();
        phiCand = candidate.phi();
        selXicCand = candidate.isSelXicToXiPiPi() >= cfgXicCand.selectionFlagXic;
        const auto& probs = candidate.mlProbXicToXiPiPi();
        fillMlOutput(probs, outputMlXic);
      }

      absPtCand = std::abs(ptCand);
      registry.fill(HIST("hPtCandXic"), absPtCand);

      phiCand = RecoDecay::constrainAngle(phiCand, -PIHalf);

      if ((std::abs(yCand) > cfgXicCand.yCandMax) || absPtCand < cfgXicCand.ptCandMin || absPtCand > cfgXicCand.ptCandMax) {
        continue;
      }

      registry.fill(HIST("hY"), yCand);
      registry.fill(HIST("hEta"), etaCand);
      registry.fill(HIST("hPhi"), phiCand);
      registry.fill(HIST("hCandBin"), poolBin);

      if (cfgXicCand.applyEfficiency) {
        efficiencyWeightCand = 1. / cfgXicCand.efficiencyXic->at(o2::analysis::findBin(cfgXicCand.binsPtEfficiencyXic, absPtCand));
      }

      if constexpr (IsMcRec) {
        isPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
        isNonPrompt = candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;

        isSignal = isXicPlus ? (std::abs(candidate.flagMcMatchRec()) == o2::aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi) : (std::abs(candidate.flagMcMatchRec()) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi)));

        if (isSignal) {
          registry.fill(HIST("hPtCandSig"), absPtCand);
          registry.fill(HIST("hEtaSig"), etaCand);
          registry.fill(HIST("hPhiSig"), phiCand);
          registry.fill(HIST("hYSig"), yCand);
        }
      }

      if (selXicCand) {
        if (isXicPlus) {
          registry.fill(HIST("hMassXicPlusVsPt"), massCand, absPtCand, efficiencyWeightCand);
          registry.fill(HIST("hMassXicPlusData"), massCand, efficiencyWeightCand);
        } else { // Xic0
          registry.fill(HIST("hMassXic0VsPt"), massCand, absPtCand, efficiencyWeightCand);
          registry.fill(HIST("hMassXic0Data"), massCand, efficiencyWeightCand);
        }

        if (isPrompt) {
          registry.fill(HIST("hPtCandSigPrompt"), absPtCand);
          registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), absPtCand, multiplicityFT0M);
        } else if (isNonPrompt) {
          registry.fill(HIST("hPtCandSigNonPrompt"), absPtCand);
          registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), absPtCand, multiplicityFT0M);
        }

        entryCandRecoInfo(massCand, ptCand, outputMlXic[0], outputMlXic[1], poolBin);
        entryCandCandGenInfo(isPrompt);
        if (!skipMixedEventTableFilling) {
          entryCand(phiCand, etaCand, ptCand, massCand, poolBin, gCollisionId, timeStamp);
          // entryCandCharge(xicType);
        }
      }

      registry.fill(HIST("hCandBin"), poolBin);

      // Correlate Xic with all Lambda V0 in the same event
      for (const auto& v0 : v0s) {
        auto posTrackV0 = v0.template posTrack_as<TrackType>();
        auto negTrackV0 = v0.template negTrack_as<TrackType>();

        if (cfgV0.cfgIsCorrCollMatchV0 && ((v0.collisionId() != posTrackV0.collisionId()) || (v0.collisionId() != negTrackV0.collisionId()))) {
          continue;
        }

        // Process Lambda (proton-pion)
        if (std::abs(o2::constants::physics::MassLambda - v0.mLambda()) < cfgV0.cfgHypMassWindow) {
          if (isSelectedV0Daughter(posTrackV0, kProton) && isSelectedV0Daughter(negTrackV0, kPiPlus)) {

            if (selXicCand) {
              fillCorrelationTable<IsMcRec, static_cast<int>(V0LambdaType::Lambda)>(cfgXicCand.fillTrkPID, v0, ptCand, etaCand, phiCand, outputMlXic, poolBin, correlationStatus, yCand, massCand, *mcParticles);
            }

            if (countCand == 1) {
              if (!skipMixedEventTableFilling) {
                entryHadron(v0.phi(), v0.eta(), v0.pt() * static_cast<int>(V0LambdaType::Lambda), poolBin, gCollisionId, timeStamp);
                entryV0InvMass(v0.mLambda());
                registry.fill(HIST("hTracksBin"), poolBin);
              }
            }
          }
        }

        // Process anti-Lambda (anti-proton-pion)
        if (std::abs(o2::constants::physics::MassLambda - v0.mAntiLambda()) < cfgV0.cfgHypMassWindow) {
          if (isSelectedV0Daughter(negTrackV0, kProton) && isSelectedV0Daughter(posTrackV0, kPiPlus)) {

            if (selXicCand) {
              fillCorrelationTable<IsMcRec, static_cast<int>(V0LambdaType::AntiLambda)>(cfgXicCand.fillTrkPID, v0, ptCand, etaCand, phiCand, outputMlXic, poolBin, correlationStatus, yCand, massCand, *mcParticles);
            }

            if (countCand == 1) {
              if (!skipMixedEventTableFilling) {
                entryHadron(v0.phi(), v0.eta(), v0.pt() * static_cast<int>(V0LambdaType::AntiLambda), poolBin, gCollisionId, timeStamp);
                entryV0InvMass(v0.mAntiLambda());
                registry.fill(HIST("hTracksBin"), poolBin);
              }
            }
          }
        }

      } // end v0 loop

      countCand++;
    } // end outer Xic loop
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), multiplicityFT0M);
  }

  template <typename T1, typename T2, typename McPart>
  void calculateTrkEff(T1 const& trackPos1, T2 const& trackPos2, McPart const& mcParticles)
  {
    decltype(trackPos1.template mcParticle_as<aod::McParticles>()) mctrk{};
    if (trackPos1.has_mcParticle()) {
      mctrk = trackPos1.template mcParticle_as<aod::McParticles>();
    } else if (trackPos2.has_mcParticle()) {
      mctrk = trackPos2.template mcParticle_as<aod::McParticles>();
    } else {
      return;
    }
    auto gentracks = mcParticles.sliceBy(perTrueCollision, mctrk.mcCollisionId());
    for (const auto& track : gentracks) {
      if (std::abs(track.eta()) > cfgXicCand.etaTrackMax || track.pt() < cfgXicCand.ptTrackMin || track.pt() > cfgXicCand.ptTrackMax) {
        continue;
      }
      if ((std::abs(track.pdgCode()) != kElectron) && (std::abs(track.pdgCode()) != kMuonMinus) && (std::abs(track.pdgCode()) != kPiPlus) && (std::abs(track.pdgCode()) != kKPlus) && (std::abs(track.pdgCode()) != kProton)) {
        continue;
      }
      if (cfgXicCand.pidTrkApplied && (std::abs(track.pdgCode()) != cfgXicCand.particlePdg)) {
        continue;
      }

      if (!track.isPhysicalPrimary()) {
        continue;
      }

      auto motherTrkGen = mcParticles.iteratorAt(track.mothersIds()[0]);
      if (std::abs(motherTrkGen.pdgCode()) == kXiC0 || std::abs(motherTrkGen.pdgCode()) == kXiCPlus) {
        continue;
      }

      auto chargeTrack = pdg->GetParticle(track.pdgCode())->Charge();
      registry.fill(HIST("hPtTracksVsSignGen"), track.pt(), chargeTrack / (std::abs(chargeTrack)));
    }
  }

  template <bool IsMcRec, int LambdaType, typename AssocType, typename McPart>
  void fillCorrelationTable(bool trkPidFill, AssocType const& assoc, float const& candPt, float const& candEta, float const& candPhi,
                            const std::vector<float>& outMl, int binPool, int8_t correlStatus,
                            double yCand, double candMass, McPart const& mcParticles)
  {
    bool isPhysicalPrimary = false;
    int trackOrigin = -1;
    float const cent = 100.f;
    double massCandHadron = -999.0;
    double ptCombined = -999.0;
    double yAssoc = -999.0;
    int signAssoc = 0;

    if constexpr (LambdaType == static_cast<int>(V0LambdaType::Lambda)) { // Lambda
      massCandHadron = calculateInvMass(candPt, candEta, candPhi, assoc, candMass, assoc.mLambda());
      entryPairedV0InvMass(assoc.mLambda());
      signAssoc = static_cast<int>(V0LambdaType::Lambda);
      yAssoc = assoc.yLambda();
    } else if constexpr (LambdaType == static_cast<int>(V0LambdaType::AntiLambda)) { // AntiLambda
      massCandHadron = calculateInvMass(candPt, candEta, candPhi, assoc, candMass, assoc.mAntiLambda());
      entryPairedV0InvMass(assoc.mAntiLambda());
      signAssoc = static_cast<int>(V0LambdaType::AntiLambda);
      yAssoc = assoc.yLambda();
    } else { // Regular track
      massCandHadron = calculateInvMass(candPt, candEta, candPhi, assoc, candMass, getMassFromPdg(cfgXicCand.particlePdg));
      signAssoc = assoc.sign();
      yAssoc = assoc.rapidity(getMassFromPdg(cfgXicCand.particlePdg));
    }

    ptCombined = calculateCombinedPt(candPt, candEta, candPhi, assoc);

    registry.fill(HIST("hXicHadronInvMassVsPt"), massCandHadron, ptCombined);
    registry.fill(HIST("hXicHadronInvMass"), massCandHadron);
    registry.fill(HIST("hXicHadronPtCombined"), ptCombined);

    entryCandHadronPair(getDeltaPhi(assoc.phi(), candPhi),
                        assoc.eta() - candEta,
                        candPt,
                        assoc.pt() * signAssoc,
                        binPool,
                        correlStatus,
                        cent); // will be added later if required
    entryCandHadronPairY(yAssoc - yCand);
    entryCandHadronMlInfo(outMl[0], outMl[1]);
    entryCandHadronRecoInfo(candMass, false);

    if constexpr (LambdaType == 0) { // Regular track
      entryTrackRecoInfo(assoc.dcaXY(), assoc.dcaZ(), assoc.tpcNClsCrossedRows());
      if (trkPidFill) {
        entryCandHadronPairTrkPID(assoc.tpcNSigmaPr(), assoc.tpcNSigmaKa(), assoc.tpcNSigmaPi(), assoc.tofNSigmaPr(), assoc.tofNSigmaKa(), assoc.tofNSigmaPi());
      }
    }
    entryXicHadronInvMass(massCandHadron, ptCombined);

    if constexpr (IsMcRec) {
      if (assoc.has_mcParticle()) {
        auto mcParticle = assoc.template mcParticle_as<aod::McParticles>();
        isPhysicalPrimary = mcParticle.isPhysicalPrimary();
        trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, true);
        entryCandHadronGenInfo(isPrompt, isPhysicalPrimary, trackOrigin);
      } else {
        entryCandHadronGenInfo(isPrompt, false, 0);
        registry.fill(HIST("hFakeTracksMcRec"), assoc.pt());
      }
      registry.fill(HIST("hPtParticleAssocVsCandMcRec"), assoc.pt(), candPt);
      if (isPhysicalPrimary) {
        registry.fill(HIST("hPtPrimaryParticleAssocVsCandMcRec"), assoc.pt(), candPt);
      }
    }
  }

  // ============================================================================
  // SAME EVENT PROCESSING (WITH REGULAR HADRON TRACKS)
  // ============================================================================

  template <bool IsMcRec, bool isXicPlus, typename CollisionType, typename CandType, typename TrackType>
  void doSameEvent(CollisionType const& collision,
                   TrackType const& tracks,
                   CandType const& candidates,
                   aod::McParticles const* mcParticles = nullptr)
  {
    int nTracks = 0;
    int64_t timeStamp = 0;
    bool skipMixedEventTableFilling = false;
    float const multiplicityFT0M = collision.multFT0M();
    int gCollisionId = collision.globalIndex();

    if (cfgXicCand.eventFractionToAnalyze > 0) {
      if (rnd->Uniform(0, 1) > cfgXicCand.eventFractionToAnalyze) {
        skipMixedEventTableFilling = true;
      }
    }

    if constexpr (!IsMcRec) {
      timeStamp = collision.template bc_as<aod::BCsWithTimestamps>().timestamp();
    }

    poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), multiplicityFT0M));

    if (cfgXicCand.correlateXicWithLeadingParticle) {
      leadingIndex = findLeadingParticle(tracks, cfgXicCand.etaTrackMax.value);
    }

    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (std::abs(track.eta()) > cfgXicCand.etaTrackMax || std::abs(track.dcaXY()) > cfgXicCand.dcaXYTrackMax || std::abs(track.dcaZ()) > cfgXicCand.dcaZTrackMax) {
          continue;
        }
        nTracks++;
      }
    }

    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < cfgXicCand.multMin || nTracks > cfgXicCand.multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    int countCand = 1;

    for (const auto& candidate : candidates) {
      double efficiencyWeightCand = 1.;
      double yCand = -999.0;
      double etaCand = -999.0;
      double ptCand = -999.0;
      double phiCand = -999.0;
      double massCand = -999.0;
      bool selXicCand = false;
      bool isCandidateDaughter = true;

      if constexpr (!isXicPlus) {
        massCand = candidate.invMassCharmBaryon();
        yCand = candidate.kfRapXic(); // yCand = candidate.y(o2::constants::physics::MassXiC0);
        ptCand = -RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon()) * candidate.signDecay();
        etaCand = candidate.etaCharmBaryon();
        phiCand = RecoDecay::phi(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
        selXicCand = candidate.resultSelections() >= cfgXicCand.selectionFlagXic;
        const auto& probs = candidate.mlProbToXiPi();
        fillMlOutput(probs, outputMlXic);
      } else { // XicPlus

        massCand = candidate.invMassXicPlus();
        yCand = candidate.y(o2::constants::physics::MassXiCPlus);
        ptCand = candidate.pt() * candidate.sign();
        etaCand = candidate.eta();
        phiCand = candidate.phi();
        selXicCand = candidate.isSelXicToXiPiPi() >= cfgXicCand.selectionFlagXic;
        const auto& probs = candidate.mlProbXicToXiPiPi();
        fillMlOutput(probs, outputMlXic);
      }

      double absPtCand = std::abs(ptCand);
      registry.fill(HIST("hPtCandXic"), absPtCand);

      if ((std::abs(yCand) > cfgXicCand.yCandMax) || absPtCand < cfgXicCand.ptCandMin || absPtCand > cfgXicCand.ptCandMax) {
        continue;
      }

      registry.fill(HIST("hY"), yCand);
      registry.fill(HIST("hEta"), etaCand);
      registry.fill(HIST("hPhi"), phiCand);
      registry.fill(HIST("hCandBin"), poolBin);

      if (cfgXicCand.applyEfficiency) {
        efficiencyWeightCand = 1. / cfgXicCand.efficiencyXic->at(o2::analysis::findBin(cfgXicCand.binsPtEfficiencyXic, absPtCand));
      }

      if constexpr (IsMcRec) {
        isPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
        isNonPrompt = candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;

        isSignal = std::abs(candidate.flagMcMatchRec()) == o2::aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi;
        // auto trackProng0 = candidate.template prong0_as<TrackType>();
        // auto trackProng1 = candidate.template prong1_as<TrackType>();
        if (cfgXicCand.calTrkEff && countCand == 1 && (isSignal || !cfgXicCand.calEffEventWithCand)) {
          // calculateTrkEff can be added here if needed
        }
      }

      if (isSignal) {
        registry.fill(HIST("hPtCandSig"), absPtCand);
        registry.fill(HIST("hEtaSig"), etaCand);
        registry.fill(HIST("hPhiSig"), phiCand);
        registry.fill(HIST("hYSig"), yCand);
      }

      if (selXicCand) {
        if (isXicPlus) {
          registry.fill(HIST("hMassXicPlusVsPt"), massCand, absPtCand, efficiencyWeightCand);
          registry.fill(HIST("hMassXicPlusData"), massCand, efficiencyWeightCand);
        } else { // Xic0
          registry.fill(HIST("hMassXic0VsPt"), massCand, absPtCand, efficiencyWeightCand);
          registry.fill(HIST("hMassXic0Data"), massCand, efficiencyWeightCand);
        }

        if (isPrompt) {
          registry.fill(HIST("hPtCandSigPrompt"), absPtCand);
          registry.fill(HIST("hPtVsMultiplicityMcRecPrompt"), absPtCand, multiplicityFT0M);
        } else if (isNonPrompt) {
          registry.fill(HIST("hPtCandSigNonPrompt"), absPtCand);
          registry.fill(HIST("hPtVsMultiplicityMcRecNonPrompt"), absPtCand, multiplicityFT0M);
        }

        entryCandRecoInfo(massCand, ptCand, outputMlXic[0], outputMlXic[1], poolBin);
        entryCandCandGenInfo(isPrompt);
        if (!skipMixedEventTableFilling) {
          entryCand(phiCand, etaCand, ptCand, massCand, poolBin, gCollisionId, timeStamp);
        }
      }

      registry.fill(HIST("hCandBin"), poolBin);

      // Correlation with hadrons
      for (const auto& track : tracks) {

        isCandidateDaughter = (candidate.bachelorId() == track.globalIndex()) || (candidate.posTrackId() == track.globalIndex()) || (candidate.negTrackId() == track.globalIndex());

        if constexpr (isXicPlus) {
          isCandidateDaughter = (candidate.pi0Id() == track.globalIndex()) || (candidate.pi1Id() == track.globalIndex()) || isCandidateDaughter;
        } else {
          isCandidateDaughter = (candidate.bachelorFromCharmBaryonId() == track.globalIndex()) || isCandidateDaughter;
        }

        if (isCandidateDaughter && !cfgXicCand.storeAutoCorrelationFlag) {
          continue;
        }
        if (isCandidateDaughter) {
          correlationStatus = true;
        }
        if (!track.isGlobalTrackWoDCA()) {
          continue;
        }

        if (cfgXicCand.pidTrkApplied) {
          if (!passPIDSelection(track, cfgXicCand.trkPIDspecies, cfgXicCand.pidTPCMax, cfgXicCand.pidTOFMax, cfgXicCand.tofPIDThreshold, cfgXicCand.forceTOF)) {
            continue;
          }
        }

        if (cfgXicCand.correlateXicWithLeadingParticle) {
          if (track.globalIndex() != leadingIndex) {
            continue;
          }
        }

        if constexpr (IsMcRec) {
          if (cfgXicCand.calTrkEff && countCand == 1 && (isSignal || !cfgXicCand.calEffEventWithCand) && track.has_mcParticle()) {
            auto mcParticle = track.template mcParticle_as<aod::McParticles>();
            if (!mcParticle.isPhysicalPrimary() && cfgXicCand.isRecTrkPhyPrimary) {
              continue;
            }

            auto motherTrk = mcParticles->iteratorAt(mcParticle.mothersIds()[0]);
            if (std::abs(motherTrk.pdgCode()) == kXiCPlus || std::abs(motherTrk.pdgCode()) == kXiC0) {
              continue;
            }

            registry.fill(HIST("hPtTracksVsSignRec"), track.pt(), track.sign());
            if (std::abs(mcParticle.pdgCode()) == cfgXicCand.particlePdg) {
              registry.fill(HIST("hPtTracksVsSignRecTrue"), track.pt(), track.sign());
            }
          } else {
            registry.fill(HIST("hFakeTracksMcRec"), track.pt());
          }
        }

        if (selXicCand) {
          fillCorrelationTable<IsMcRec, 0>(cfgXicCand.fillTrkPID, track, ptCand, etaCand, phiCand, outputMlXic, poolBin, correlationStatus, yCand, massCand, *mcParticles);
        }

        if (countCand == 1) {
          if (!skipMixedEventTableFilling) {
            entryHadron(track.phi(), track.eta(), track.pt() * track.sign(), poolBin, gCollisionId, timeStamp);
            if (cfgXicCand.fillTrkPID) {
              entryTrkPID(track.tpcNSigmaPr(), track.tpcNSigmaKa(), track.tpcNSigmaPi(), track.tofNSigmaPr(), track.tofNSigmaKa(), track.tofNSigmaPi());
            }
            registry.fill(HIST("hTracksBin"), poolBin);
          }
        }
      } // end Hadron Tracks loop
      countCand++;
    } // end outer Xic loop
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hMultFT0M"), multiplicityFT0M);
  }
  // ============================================================================
  // MIXED EVENT PROCESSING
  // ============================================================================

  template <bool IsMcRec, bool isXicPlus, bool isV0, typename CollisionType, typename CandType, typename AssociateType, typename TrackType>
  void doMixEvent(CollisionType const& collisions,
                  AssociateType const& tracks,
                  CandType const& candidates,
                  TrackType const&,
                  aod::McParticles const* mcParticles = nullptr)
  {
    if (candidates.size() == 0) {
      return;
    }

    double yCand = -999.;
    double ptCand = -999.;
    double phiCand = -999.;
    double etaCand = -999.;
    double massCand = -999.0;
    bool selCand = false;

    auto tracksTuple = std::make_tuple(candidates, tracks);
    Pair<CollisionType, CandType, AssociateType, BinningType> const pairData{corrBinning, cfgXicCand.numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairData) {
      poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
      poolBinXic = corrBinning.getBin(std::make_tuple(c1.posZ(), c1.multFT0M()));
      registry.fill(HIST("hMultFT0M"), c1.multFT0M());
      registry.fill(HIST("hZvtx"), c1.posZ());
      registry.fill(HIST("hTracksPoolBin"), poolBin);
      registry.fill(HIST("hXicPoolBin"), poolBinXic);

      for (const auto& [candidate, assocParticle] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if constexpr (!isXicPlus) {
          massCand = candidate.invMassCharmBaryon();
          yCand = candidate.kfRapXic();
          ptCand = -RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon()) * candidate.signDecay();
          etaCand = candidate.etaCharmBaryon();
          phiCand = RecoDecay::phi(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
          selCand = candidate.resultSelections() >= cfgXicCand.selectionFlagXic;
          const auto& probs = candidate.mlProbToXiPi();
          fillMlOutput(probs, outputMlXic);
        } else {
          massCand = candidate.invMassXicPlus();
          yCand = candidate.y(o2::constants::physics::MassXiCPlus);
          ptCand = candidate.pt() * candidate.sign();
          etaCand = candidate.eta();
          phiCand = candidate.phi();
          selCand = candidate.isSelXicToXiPiPi() >= cfgXicCand.selectionFlagXic;
          const auto& probs = candidate.mlProbXicToXiPiPi();
          fillMlOutput(probs, outputMlXic);
        }

        if (!selCand) {
          continue;
        }
        if constexpr (IsMcRec) {
          isPrompt = candidate.originMcRec() == RecoDecay::OriginType::Prompt;
          isNonPrompt = candidate.originMcRec() == RecoDecay::OriginType::NonPrompt;
        }

        if (std::abs(yCand) > cfgXicCand.yCandMax) {
          continue;
        }

        if constexpr (!isV0) {
          if (cfgXicCand.corrParticle != 2 && !assocParticle.isGlobalTrackWoDCA()) {
            continue;
          }
          if (cfgXicCand.pidTrkApplied) {
            if (!passPIDSelection(assocParticle, cfgXicCand.trkPIDspecies, cfgXicCand.pidTPCMax, cfgXicCand.pidTOFMax, cfgXicCand.tofPIDThreshold, cfgXicCand.forceTOF)) {
              continue;
            }
          }
        } else {

          auto posTrackV0 = assocParticle.template posTrack_as<TrackType>();
          auto negTrackV0 = assocParticle.template negTrack_as<TrackType>();

          if (std::abs(o2::constants::physics::MassLambda - assocParticle.mLambda()) < cfgV0.cfgHypMassWindow) {
            if (isSelectedV0Daughter(posTrackV0, kProton) && isSelectedV0Daughter(negTrackV0, kPiPlus)) {

              fillCorrelationTable<IsMcRec, static_cast<int>(V0LambdaType::Lambda)>(cfgXicCand.fillTrkPID, assocParticle, ptCand, etaCand, phiCand, outputMlXic, poolBin, correlationStatus, yCand, massCand, *mcParticles);
            }
          }

          if (std::abs(o2::constants::physics::MassLambda - assocParticle.mAntiLambda()) < cfgV0.cfgHypMassWindow) {
            if (isSelectedV0Daughter(negTrackV0, kProton) && isSelectedV0Daughter(posTrackV0, kPiPlus)) {
              fillCorrelationTable<IsMcRec, static_cast<int>(V0LambdaType::AntiLambda)>(cfgXicCand.fillTrkPID, assocParticle, ptCand, etaCand, phiCand, outputMlXic, poolBin, correlationStatus, yCand, massCand, *mcParticles);
            }
          }
        }
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

  // ============================================================================
  // MC GEN LEVEL PROCESSING
  // ============================================================================

  template <bool isXicPlus, typename CollisionType, typename PartType>
  void doSameEventMcGen(CollisionType const& mcCollision, PartType const& mcParticles)
  {
    int counterCharmCand = 0;
    int8_t candSign = 0;

    registry.fill(HIST("hMcEvtCount"), 0);
    BinningTypeMcGen const corrBinningMcGen{{binsZVtx, binsMultiplicityMc}, true};
    poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), mcCollision.multMCFT0A()));
    registry.fill(HIST("hMultFT0AMcGen"), mcCollision.multMCFT0A());

    // find leading particle
    if (cfgXicCand.correlateXicWithLeadingParticle) {
      leadingIndex = findLeadingParticleMcGen(mcParticles, cfgXicCand.etaTrackMax.value, cfgXicCand.ptTrackMin.value);
    }

    // Mc Gen level
    for (const auto& particle : mcParticles) {
      if ((isXicPlus && std::abs(particle.pdgCode()) != kXiCPlus) || (!isXicPlus && std::abs(particle.pdgCode()) != kXiC0)) {
        continue;
      }

      if (particle.pdgCode() > 0) {
        candSign = static_cast<int8_t>(ParticleType::Particle);
      } else {
        candSign = static_cast<int8_t>(ParticleType::AntiParticle);
      }

      double const massCand = (std::abs(particle.pdgCode()) == kXiC0) ? o2::constants::physics::MassXiC0 : o2::constants::physics::MassXiCPlus;
      double const yCand = RecoDecay::y(particle.pVector(), massCand);

      if (std::abs(yCand) > cfgXicCand.yCandGenMax || particle.pt() < cfgXicCand.ptCandMin) {
        continue;
      }

      registry.fill(HIST("hCandBin"), poolBin);
      registry.fill(HIST("hPtCandMcGen"), particle.pt());
      registry.fill(HIST("hEtaMcGen"), particle.eta());
      registry.fill(HIST("hPhiMcGen"), RecoDecay::constrainAngle(particle.phi(), -PIHalf));
      registry.fill(HIST("hYMcGen"), yCand);

      // int8_t chargeCand = pdg->GetParticle(particle.pdgCode())->Charge() / PDGChargeScale;
      float xicType = getXicTypeMC(particle);

      isPrompt = particle.originMcGen() == RecoDecay::OriginType::Prompt;
      isNonPrompt = particle.originMcGen() == RecoDecay::OriginType::NonPrompt;

      if (isPrompt) {
        registry.fill(HIST("hPtCandMcGenPrompt"), particle.pt());
        registry.fill(HIST("hPtCandVsChargeMcGenPrompt"), particle.pt(), xicType);
      } else if (isNonPrompt) {
        registry.fill(HIST("hPtCandMcGenNonPrompt"), particle.pt());
        registry.fill(HIST("hPtCandVsChargeMcGenNonPrompt"), particle.pt(), xicType);
      }

      // Xic Hadron correlation dedicated section
      registry.fill(HIST("hcountCandtriggersMcGen"), 0, particle.pt());

      static constexpr std::size_t NDaughtersXic0 = static_cast<std::size_t>(XicDecayDaughtersCount::Xic0DaughtersCount);
      static constexpr std::size_t NDaughtersXicPlus = static_cast<std::size_t>(XicDecayDaughtersCount::XicPlusDaughtersCount);
      std::vector<int> listDaughters{};
      listDaughters.clear();
      const std::size_t nDaughtersExpected = (std::abs(xicType) > 1.) ? NDaughtersXicPlus : NDaughtersXic0;

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

      for (const auto& particleAssoc : mcParticles) {
        if (std::abs(particleAssoc.eta()) > cfgXicCand.etaTrackMax || particleAssoc.pt() < cfgXicCand.ptTrackMin || particleAssoc.pt() > cfgXicCand.ptTrackMax) {
          continue;
        }

        if (cfgXicCand.corrParticle == 1 && ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particleAssoc.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton))) {
          continue;
        }

        if (cfgXicCand.corrParticle == 2 && (std::abs(particleAssoc.pdgCode()) != cfgXicCand.particlePdg)) {
          continue;
        }

        if (cfgXicCand.corrParticle == 1 && !particleAssoc.isPhysicalPrimary()) {
          continue;
        }

        if (cfgXicCand.correlateXicWithLeadingParticle) {
          if (particleAssoc.globalIndex() != leadingIndex) {
            continue;
          }
        }

        if (std::find(prongsId.begin(), prongsId.end(), particleAssoc.globalIndex()) != prongsId.end()) {
          if (!cfgXicCand.storeAutoCorrelationFlag) {
            continue;
          }
          correlationStatus = true;
        }

        int trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, particleAssoc, true);
        int8_t chargeAssoc = pdg->GetParticle(particleAssoc.pdgCode())->Charge();
        chargeAssoc = chargeAssoc / std::abs(chargeAssoc);
        registry.fill(HIST("hPtParticleAssocMcGen"), particleAssoc.pt());
        float const cent = 100.0;

        entryCandHadronPair(getDeltaPhi(particleAssoc.phi(), particle.phi()),
                            particleAssoc.eta() - particle.eta(),
                            particle.pt() * candSign,
                            particleAssoc.pt() * chargeAssoc,
                            poolBin,
                            correlationStatus,
                            cent);
        entryCandHadronPairY(RecoDecay::y(particleAssoc.pVector(), getMassFromPdg(particleAssoc.pdgCode())) - yCand);
        entryCandHadronRecoInfo(massCand, true);
        entryCandHadronGenInfo(isPrompt, particleAssoc.isPhysicalPrimary(), trackOrigin);
      } // end inner loop
    } // end outer loop
    registry.fill(HIST("hcountCandHadronPerEvent"), counterCharmCand);
    registry.fill(HIST("hZvtx"), mcCollision.posZ());
  }

  // ============================================================================
  // PROCESS FUNCTIONS
  // ============================================================================

  /// Data processing: XicPlus with regular hadron tracks
  void processDataXicPlus(SelCollisions::iterator const& collision,
                          TracksData const& tracks,
                          CandsXicPlusDataFiltered const& candidates,
                          aod::BCsWithTimestamps const&)
  {

    doSameEvent<false, 1>(collision, tracks, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processDataXicPlus, "Process data XicPlus", true);

  /// MC Reco processing: XicPlus with regular hadron tracks
  void processMcRecXicPlus(SelCollisions::iterator const& collision,
                           TracksWithMc const& tracks,
                           CandsXicPlusMcRecFiltered const& candidates,
                           aod::McParticles const& mcParticles)
  {
    doSameEvent<true, 1>(collision, tracks, candidates, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcRecXicPlus, "Process Mc Reco mode Xic", false);

  /// Data processing: Xic0 with regular hadron tracks
  void processDataXic0(SelCollisions::iterator const& collision,
                       TracksData const& tracks,
                       CandsXic0DataFiltered const& candidates,
                       aod::BCsWithTimestamps const&)
  {
    doSameEvent<false, 0>(collision, tracks, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processDataXic0, "Process data Xic0", false);

  /// MC Reco processing: Xic0 with regular hadron tracks
  void processMcRecXic0(SelCollisions::iterator const& collision,
                        TracksWithMc const& tracks,
                        CandsXic0McRecFiltered const& candidates,
                        aod::McParticles const& mcParticles)
  {
    doSameEvent<true, 0>(collision, tracks, candidates, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcRecXic0, "Process Mc Reco mode Xic0", false);

  /// Data processing: XicPlus with V0 Lambda
  void processDataXicPlusV0(SelCollisions::iterator const& collision,
                            TracksData const& tracks,
                            aod::V0Datas const& v0s,
                            CandsXicPlusDataFiltered const& candidates,
                            aod::BCsWithTimestamps const&)
  {
    doSameEventWithV0<false, 1>(collision, v0s, tracks, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processDataXicPlusV0, "Data process for v0 lambda with Xic Plus", false);

  /// MC Reco processing: XicPlus with V0 Lambda
  void processMcRecXicPlusV0(SelCollisions::iterator const& collision,
                             TracksWithMc const& tracks,
                             soa::Join<aod::V0Datas, aod::McV0Labels> const& v0s,
                             CandsXicPlusMcRecFiltered const& candidates,
                             aod::McParticles const& mcParticles)
  {
    doSameEventWithV0<true, 1>(collision, v0s, tracks, candidates, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcRecXicPlusV0, "Mc process for v0 lambda with Xic Plus", false);

  /// Data processing: Xic0 with V0 Lambda
  void processDataXic0V0(SelCollisions::iterator const& collision,
                         TracksData const& tracks,
                         aod::V0Datas const& v0s,
                         CandsXic0DataFiltered const& candidates,
                         aod::BCsWithTimestamps const&)
  {
    doSameEventWithV0<false, 0>(collision, v0s, tracks, candidates);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processDataXic0V0, "Data process for v0 lambda with Xic0", false);

  /// MC Reco processing: Xic0 with V0 Lambda
  void processMcRecXic0V0(SelCollisions::iterator const& collision,
                          TracksWithMc const& tracks,
                          soa::Join<aod::V0Datas, aod::McV0Labels> const& v0s,
                          CandsXic0McRecFiltered const& candidates,
                          aod::McParticles const& mcParticles)
  {
    doSameEventWithV0<true, 0>(collision, v0s, tracks, candidates, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcRecXic0V0, "Mc process for v0 lambda with Xic0", false);

  // ============================================================================
  // MIXED EVENT PROCESS FUNCTIONS
  // ============================================================================

  /// Data Mixed Event: XicPlus with regular hadron tracks
  void processDataMixedEventXicPlus(SelCollisions const& collisions,
                                    CandsXicPlusDataFiltered const& candidates,
                                    TracksData const& tracks)
  {
    doMixEvent<false, 1, 0>(collisions, tracks, candidates, tracks);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processDataMixedEventXicPlus, "Process Mixed Event Data XicPlus", false);

  /// MC Reco Mixed Event: XicPlus with regular hadron tracks
  void processMcRecMixedEventXicPlus(SelCollisions const& collisions,
                                     CandsXicPlusMcRecFiltered const& candidates,
                                     TracksWithMc const& tracks,
                                     aod::McParticles const& mcParticles)
  {
    doMixEvent<true, 1, 0>(collisions, tracks, candidates, tracks, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcRecMixedEventXicPlus, "Process Mixed Event McRec XicPlus", false);

  /// Data Mixed Event: Xic0 with regular hadron tracks
  void processDataMixedEventXic0(SelCollisions const& collisions,
                                 CandsXic0DataFiltered const& candidates,
                                 TracksData const& tracks)
  {
    doMixEvent<false, 0, 0>(collisions, tracks, candidates, tracks);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processDataMixedEventXic0, "Process Mixed Event Data Xic0", false);

  /// MC Reco Mixed Event: Xic0 with regular hadron tracks
  void processMcRecMixedEventXic0(SelCollisions const& collisions,
                                  CandsXic0McRecFiltered const& candidates,
                                  TracksWithMc const& tracks,
                                  aod::McParticles const& mcParticles)
  {
    doMixEvent<true, 0, 0>(collisions, tracks, candidates, tracks, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcRecMixedEventXic0, "Process Mixed Event McRec Xic0", false);

  // ============================================================================
  // MIXED EVENT WITH V0 LAMBDA - NEW FUNCTIONS
  // ============================================================================

  /// Data Mixed Event: XicPlus with V0 Lambda
  /// NOTE: V0 mixed events are more complex - need proper binning and collision matching
  void processDataMixedEventXicPlusV0(SelCollisions const& collisions,
                                      CandsXicPlusDataFiltered const& candidates,
                                      aod::V0Datas const& v0s,
                                      TracksData const& tracks)
  {
    doMixEvent<false, 1, 1>(collisions, v0s, candidates, tracks);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processDataMixedEventXicPlusV0, "Process Mixed Event Data XicPlus + V0Lambda", false);

  /// MC Reco Mixed Event: XicPlus with V0 Lambda
  void processMcRecMixedEventXicPlusV0(SelCollisions const& collisions,
                                       CandsXicPlusMcRecFiltered const& candidates,
                                       soa::Join<aod::V0Datas, aod::McV0Labels> const& v0s,
                                       TracksWithMc const& tracks,
                                       aod::McParticles const& mcParticles)
  {
    doMixEvent<true, 1, 1>(collisions, v0s, candidates, tracks, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcRecMixedEventXicPlusV0, "Process Mixed Event McRec XicPlus + V0Lambda", false);

  /// Data Mixed Event: Xic0 with V0 Lambda
  void processDataMixedEventXic0V0(SelCollisions const& collisions,
                                   CandsXic0DataFiltered const& candidates,
                                   aod::V0Datas const& v0s,
                                   TracksData const& tracks)
  {
    doMixEvent<false, 0, 1>(collisions, v0s, candidates, tracks);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processDataMixedEventXic0V0, "Process Mixed Event Data Xic0 + V0Lambda", false);

  /// MC Reco Mixed Event: Xic0 with V0 Lambda
  void processMcRecMixedEventXic0V0(SelCollisions const& collisions,
                                    CandsXic0McRecFiltered const& candidates,
                                    soa::Join<aod::V0Datas, aod::McV0Labels> const& v0s,
                                    TracksWithMc const& tracks,
                                    aod::McParticles const& mcParticles)
  {
    doMixEvent<true, 0, 1>(collisions, v0s, candidates, tracks, &mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcRecMixedEventXic0V0, "Process Mixed Event McRec Xic0 + V0Lambda", false);

  // ============================================================================
  // MC GEN LEVEL - SAME EVENT PROCESSING
  // ============================================================================

  /// MC Gen Same Event: XicPlus
  void processMcGenXicPlus(SelCollisionsMc::iterator const& mcCollision,
                           CandidatesXicPlusMcGen const& mcParticles)
  {
    doSameEventMcGen<1>(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcGenXicPlus, "Process Mc Gen XicPlus", false);

  /// MC Gen Same Event: Xic0
  void processMcGenXic0(SelCollisionsMc::iterator const& mcCollision,
                        CandidatesXic0McGen const& mcParticles)
  {
    doSameEventMcGen<0>(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcGenXic0, "Process Mc Gen Xic0", false);

  void processMcGenMixedEvent(SelCollisionsMc const& collisions,
                              CandidatesXicPlusMcGen const& mcParticles)
  {
    BinningTypeMcGen const corrBinningMcGen{{binsZVtx, binsMultiplicityMc}, true};
    auto tracksTuple = std::make_tuple(mcParticles, mcParticles);
    Pair<SelCollisionsMc, CandidatesXicPlusMcGen, CandidatesXicPlusMcGen, BinningTypeMcGen> const pairMcGen{corrBinningMcGen, cfgXicCand.numberEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairMcGen) {
      poolBin = corrBinningMcGen.getBin(std::make_tuple(c1.posZ(), c1.multMCFT0A()));

      for (const auto& [candidate, particleAssoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (std::abs(candidate.pdgCode()) != kXiCPlus) {
          continue;
        }

        double const massCand = (std::abs(candidate.pdgCode()) == kXiC0) ? o2::constants::physics::MassXiC0 : o2::constants::physics::MassXiCPlus;
        double const yXic = RecoDecay::y(candidate.pVector(), massCand);
        if (std::abs(yXic) > cfgXicCand.yCandGenMax || candidate.pt() < cfgXicCand.ptCandMin || candidate.pt() > cfgXicCand.ptCandMax) {
          continue;
        }

        if (std::abs(particleAssoc.eta()) > cfgXicCand.etaTrackMax || particleAssoc.pt() < cfgXicCand.ptTrackMin || particleAssoc.pt() > cfgXicCand.ptTrackMax) {
          continue;
        }

        if ((std::abs(particleAssoc.pdgCode()) != kElectron) && (std::abs(particleAssoc.pdgCode()) != kMuonMinus) && (std::abs(particleAssoc.pdgCode()) != kPiPlus) && (std::abs(particleAssoc.pdgCode()) != kKPlus) && (std::abs(particleAssoc.pdgCode()) != kProton)) {
          continue;
        }

        if (!particleAssoc.isPhysicalPrimary()) {
          continue;
        }

        if (cfgXicCand.pidTrkApplied && (std::abs(particleAssoc.pdgCode()) != cfgXicCand.particlePdg)) {
          continue;
        }

        // int8_t const chargeXic = pdg->GetParticle(candidate.pdgCode())->Charge();
        float xicType = getXicTypeMC(candidate);

        int8_t const chargeAssoc = pdg->GetParticle(particleAssoc.pdgCode())->Charge();
        float cent = 100.0;

        int trackOrigin = RecoDecay::getCharmHadronOrigin(mcParticles, particleAssoc, true);
        bool isPromptXic = candidate.originMcGen() == RecoDecay::OriginType::Prompt;

        entryCandHadronPair(getDeltaPhi(particleAssoc.phi(), candidate.phi()),
                            particleAssoc.eta() - candidate.eta(),
                            candidate.pt(),
                            particleAssoc.pt() * chargeAssoc / std::abs(chargeAssoc),
                            poolBin,
                            correlationStatus,
                            cent);
        entryCandHadronPairY(RecoDecay::y(particleAssoc.pVector(), getMassFromPdg(particleAssoc.pdgCode())) - yXic); // particleAssoc.rapidity(getMassFromPdg(particleAssoc.pdgCode()))
        entryCandHadronRecoInfo(massCand, true);
        entryCandHadronGenInfo(isPromptXic, particleAssoc.isPhysicalPrimary(), trackOrigin);
        entryPairCandCharge(xicType);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorXicHadrons, processMcGenMixedEvent, "Process Mixed Event McGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorXicHadronsSelection>(cfgc),
                      adaptAnalysisTask<HfCorrelatorXicHadrons>(cfgc)};
}
