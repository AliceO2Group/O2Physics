// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamPairTaskCharmHadrons.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Ravindra SIngh, GSI, ravindra.singh@cern.ch

#include <vector>

#include "Framework/Expressions.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/Core/femtoDreamContainer.h"
#include "PWGCF/FemtoDream/Core/femtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct femtoDreamPairTaskCharmHadrons {

  enum PairCharge {
    positiveCharge = 1,
    negativeCharge = -1
  };

  /// General options
  Configurable<bool> ConfOptCPRPlotPerRadii{"ConfOptCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<bool> ConfOptUseCPR{"ConfOptUseCPR", false, "Close Pair Rejection"};
  Configurable<float> ConfOptCPRdeltaEtaMax{"ConfOptCPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
  Configurable<float> ConfOptCPRdeltaPhiMax{"ConfOptCPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<bool> ConfOptExtendedPlots{"ConfOptExtendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
  Configurable<float> ConfOptHighkstarCut{"ConfOptHighkstarCut", 100000., "Set a cut for high k*, above which the pairs are rejected"};
  Configurable<bool> ConfOptisMc{"ConfOptisMc", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<bool> ConfOptsmearingByOrigin{"ConfOptsmearingByOrigin", false, "Obtain the smearing matrix differential in the MC origin of particle 1 and particle 2. High memory consumption. Use with care!"};
  Configurable<bool> ConfOptUse4D{"ConfOptUse4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
  ConfigurableAxis ConfOptDummy{"ConfOptDummy", {1, 0, 1}, "Dummy axis"};

  /// Event selection
  Configurable<int> ConfEvent_maxMult{"ConfEvent_maxMult", 99999, "Maximum Multiplicity (MultNtr)"};
  Configurable<float> ConfEvent_maxMultPercentile{"ConfEvent_maxMultPercentile", 100, "Maximum Multiplicity Percentile"};
  Configurable<int> ConfEvent_minMult{"ConfEvent_minMult", 0, "Minimum Multiplicity (MultNtr)"};
  Configurable<float> ConfEvent_minMultPercentile{"ConfEvent_minMultPercentile", 0, "Minimum Multiplicity Percentile"};

  /// Particle 1 (track)
  Configurable<femtodreamparticle::cutContainerType> ConfTrk1_CutBit{"ConfTrk1_CutBit", 5542474, "Particle 1 (Track) - Selection bit from cutCulator"};
  Configurable<int> ConfTrk1_PDGCode{"ConfTrk1_PDGCode", 2212, "PDG code of Particle 1 (Track)"};
  Configurable<float> ConfTrk1_PIDThres{"ConfTrk1_PIDThres", 0.75, "Momentum threshold for PID selection for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk1_TPCBit{"ConfTrk1_TPCBit", 4, "PID TPC bit from cutCulator for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk1_TPCTOFBit{"ConfTrk1_TPCTOFBit", 2, "PID TPCTOF bit from cutCulator for particle 1 (Track)"};
  Configurable<float> ConfTrk1_maxEta{"ConfTrk1_maxEta", 10., "Maximum eta of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_maxPt{"ConfTrk1_maxPt", 999., "Maximum pT of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_minEta{"ConfTrk1_minEta", -10., "Minimum eta of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_minPt{"ConfTrk1_minPt", 0., "Minimum pT of partricle 1 (Track)"};

  /// Particle 2 (Charm Hadrons)
  Configurable<float> ConfHF_bkgBDT{"ConfHF_bkgBDT", 1., "Maximum background bdt score for Charm Hadron (particle 2)"};
  Configurable<int8_t> ConfHF_CandSel{"ConfHF_CandSel", 1, "candidate selection for charm hadron"};
  Configurable<float> ConfHF_fdBDT{"ConfHF_fdBDT", 0., "Minimum feed-down bdt score Charm Hadron (particle 2)"};
  Configurable<float> ConfHF_maxInvMass{"ConfHF_maxInvMass", 2.45, "Maximum invariant mass of Charm Hadron (particle 2)"};
  Configurable<float> ConfHF_maxPt{"ConfHF_maxPt", 999., "Maximum pT of Charm Hadron (particle 2)"};
  Configurable<float> ConfHF_minInvMass{"ConfHF_minInvMass", 2.15, "Minimum invariant mass of Charm Hadron (particle 2)"};
  Configurable<float> ConfHF_minPt{"ConfHF_minPt", 0., "Minimum pT of Charm Hadron (particle 2)"};
  Configurable<int> ConfHF_PDGCode{"ConfHF_PDGCode", 4122, "PDG code of particle 2 Charm Hadron"};
  Configurable<float> ConfHF_promptBDT{"ConfHF_promptBDT", 0., "Minimum prompt bdt score Charm Hadron (particle 2)"};

  /// Binning configurables
  ConfigurableAxis ConfBin4Dkstar{"ConfBin4Dkstar", {1500, 0., 6.}, "binning kstar for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
  ConfigurableAxis ConfBin4DMult{"ConfBin4Dmult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
  ConfigurableAxis ConfBin4DmT{"ConfBin4DmT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
  ConfigurableAxis ConfBin4DmultPercentile{"ConfBin4DmultPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
  ConfigurableAxis ConfBinInvMass{"ConfBinInvMass", {300, 2.15, 2.45}, "InvMass binning"};
  ConfigurableAxis ConfBinTempFitVarHF{"ConfBinTempFitVarHF", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis ConfBinTempFitVarHFChild{"ConfBinTempFitVarHFChild", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0 child)"};
  ConfigurableAxis ConfBinTempFitVarTrack{"ConfBinTempFitVarTrack", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Track)"};
  ConfigurableAxis ConfBinmT{"ConfBinmT", {225, 0., 7.5}, "binning mT"};
  ConfigurableAxis ConfBinmultTempFit{"ConfBinmultTempFit", {1, 0, 1}, "multiplicity Binning for the TempFitVar plot"};
  ConfigurableAxis ConfBinpT{"ConfBinpT", {20, 0.5, 4.05}, "pT binning"};
  ConfigurableAxis ConfBinpTHF{"ConfBinpTHF", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis ConfBinpTHFChild{"ConfBinpTHFChild", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis ConfBinpTTrack{"ConfBinpTTrack", {50, 0.5, 10.05}, "pT binning of the pT vs. TempFitVar plot (Track)"};
  ConfigurableAxis ConfBinkT{"ConfBinkT", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfBinkstar{"ConfBinkstar", {1500, 0., 6.}, "binning kstar"};

  // Mixing configurables
  ConfigurableAxis ConfMixingBinMult{"ConfMixingBinMult", {VARIABLE_WIDTH, 0.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfMixingBinMultPercentile{"ConfMixingBinMultPercentile", {VARIABLE_WIDTH, 0.0f, 100.f}, "Mixing bins - multiplicity percentile"};
  ConfigurableAxis ConfMixingBinVztx{"ConfMixingBinVztx", {VARIABLE_WIDTH, -10.0f, -4.f, 0.f, 4.f, 10.f}, "Mixing bins - z-vertex"};
  Configurable<int> ConfMixingDepth{"ConfMixingDepth", 5, "Number of events for mixing"};
  Configurable<int> ConfMixingPolicy{"ConfMixingBinPolicy", 0, "Binning policy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{ConfMixingBinVztx, ConfMixingBinMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{ConfMixingBinVztx, ConfMixingBinMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{ConfMixingBinVztx, ConfMixingBinMult, ConfMixingBinMultPercentile}, true};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron> pairCloseRejection;

  Filter trackEtaFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta < ConfTrk1_maxEta, true);
  Filter trackEtaFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta > ConfTrk1_minEta, true);
  Filter trackPtFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt < ConfTrk1_maxPt, true);
  Filter trackPtFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt > ConfTrk1_minPt, true);
  Filter hfCandSelFilter = aod::fdhf::candidateSelFlag >= ConfHF_CandSel.value;

  using FilteredCharmCands = soa::Filtered<aod::FDHfCand>;
  using FilteredCharmCand = FilteredCharmCands::iterator;

  using FilteredCharmMCCands = soa::Filtered<soa::Join<aod::FDHfCand, aod::FDHfCandMC, aod::FDHfCandMCGen>>;
  using FilteredCharmMCCand = FilteredCharmMCCands::iterator;

  using FilteredColisions = FDCollisions;
  using FilteredColision = FilteredColisions::iterator;

  using FilteredFDMCParts = soa::Filtered<soa::Join<aod::FDParticles, aod::FDParticlesIndex, aod::FDMCLabels>>;
  using FilteredFDMCPart = FilteredFDMCParts::iterator;

  using FilteredFDParticles = soa::Filtered<soa::Join<aod::FDParticles, aod::FDParticlesIndex>>;
  using FilteredFDParticle = FilteredFDParticles::iterator;

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Partition for particle 1

  Partition<FilteredFDParticles> PartitionTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack));

  Partition<FilteredFDMCParts> PartitionMCTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 (ncheckbit(aod::femtodreamparticle::cut, ConfTrk1_CutBit)) &&
                                                 ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk1_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCTOFBit));

  /// Partition for particle 2
  Partition<FilteredCharmCands> PartitionHF = aod::fdhf::bdtBkg < ConfHF_bkgBDT && aod::fdhf::bdtPrompt > ConfHF_promptBDT;
  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"CorrelationsHF", {}, OutputObjHandlingPolicy::AnalysisObject};

  float MassOne = o2::analysis::femtoDream::getMass(ConfTrk1_PDGCode);
  float MassTwo = o2::analysis::femtoDream::getMass(ConfHF_PDGCode);
  int8_t partSign = 0;
  int64_t processType = 0;

  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  Produces<o2::aod::FDResultsHF> fillFemtoResult;

  void init(InitContext& /*context*/)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry, ConfBinmultTempFit, ConfOptDummy, ConfBinpTTrack, ConfOptDummy, ConfOptDummy, ConfBinTempFitVarTrack, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptisMc, ConfTrk1_PDGCode);

    sameEventCont.init<true>(&resultRegistry,
                             ConfBinkstar, ConfBinpT, ConfBinkT, ConfBinmT, ConfMixingBinMult, ConfMixingBinMultPercentile,
                             ConfBin4Dkstar, ConfBin4DmT, ConfBin4DMult, ConfBin4DmultPercentile,
                             ConfOptisMc, ConfOptUse4D, ConfOptExtendedPlots,
                             ConfOptHighkstarCut,
                             ConfOptsmearingByOrigin, ConfBinInvMass);

    sameEventCont.setPDGCodes(ConfTrk1_PDGCode, ConfHF_PDGCode);
    mixedEventCont.init<true>(&resultRegistry,
                              ConfBinkstar, ConfBinpT, ConfBinkT, ConfBinmT, ConfMixingBinMult, ConfMixingBinMultPercentile,
                              ConfBin4Dkstar, ConfBin4DmT, ConfBin4DMult, ConfBin4DmultPercentile,
                              ConfOptisMc, ConfOptUse4D, ConfOptExtendedPlots,
                              ConfOptHighkstarCut,
                              ConfOptsmearingByOrigin, ConfBinInvMass);

    mixedEventCont.setPDGCodes(ConfTrk1_PDGCode, ConfHF_PDGCode);
    pairCleaner.init(&qaRegistry);
    if (ConfOptUseCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfOptCPRdeltaPhiMax.value, ConfOptCPRdeltaEtaMax.value, ConfOptCPRPlotPerRadii.value);
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  template <bool isMc, typename PartitionType, typename CandType, typename TableTracks, typename TableCandidates, typename Collision>
  void doSameEvent(PartitionType& SliceTrk1, CandType& SliceCharmHad, TableTracks const& parts, TableCandidates const& /*candidates*/, Collision const& col)
  {
    processType = 1; // for same event
    /// Histogramming same event
    for (auto const& part : SliceTrk1) {

      trackHistoPartOne.fillQA<isMc, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
    }
    for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceCharmHad))) {
      // proton track charge
      float chargeTrack = 0.;
      if ((p1.cut() & 1) == 1) {
        chargeTrack = positiveCharge;
      } else {
        chargeTrack = negativeCharge;
      }

      if (chargeTrack != p2.charge())
        continue;
      float kstar = FemtoDreamMath::getkstar(p1, MassOne, p2, MassTwo);
      if (kstar > ConfOptHighkstarCut)
        continue;

      if (chargeTrack == 1) {
        partSign = 1;
      } else {
        partSign = 1 << 1;
      }

      if (ConfOptUseCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, col.magField())) {
          continue;
        }
      }

      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      float invMass;
      if (p2.candidateSelFlag() == 1) {
        invMass = p2.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
      } else {
        invMass = p2.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
      }

      if (invMass < ConfHF_minInvMass || invMass > ConfHF_maxInvMass)
        continue;
      if (p2.pt() < ConfHF_minPt || p2.pt() > ConfHF_maxPt)
        continue;

      fillFemtoResult(
        invMass,
        p2.pt(),
        p1.pt(),
        p2.bdtBkg(),
        p2.bdtPrompt(),
        p2.bdtFD(),
        kstar,
        FemtoDreamMath::getkT(p1, MassOne, p2, MassTwo),
        FemtoDreamMath::getmT(p1, MassOne, p2, MassTwo),
        col.multNtr(),
        col.multV0M(),
        partSign,
        processType);

      sameEventCont.setPair<isMc, true>(p1, p2, col.multNtr(), col.multV0M(), ConfOptUse4D, ConfOptExtendedPlots, ConfOptsmearingByOrigin);
    }
  }

  void processSameEvent(FilteredColision const& col, FilteredFDParticles const& parts, FilteredCharmCands const& candidates)
  {
    eventHisto.fillQA(col);
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceCharmHad = PartitionHF->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<false>(SliceTrk1, SliceCharmHad, parts, candidates, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskCharmHadrons, processSameEvent, "Enable processing same event", true);

  template <bool isMc, typename CollisionType, typename PartType, typename PartitionType1, typename PartitionType2, typename BinningType>
  void doMixedEvent(CollisionType const& cols, PartType const& parts, PartitionType1& part1, PartitionType2& /*part2*/, BinningType policy)
  {
    processType = 1 << 1; // for mixed event

    for (auto const& [collision1, collision2] : soa::selfCombinations(policy, ConfMixingDepth.value, -1, cols, cols)) {

      auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceCharmHad = PartitionHF->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceCharmHad))) {

        float chargeTrack = 0.;
        if ((p1.cut() & 1) == 1) {
          chargeTrack = positiveCharge;
        } else {
          chargeTrack = negativeCharge;
        }

        if (chargeTrack != p2.charge())
          continue;
        float kstar = FemtoDreamMath::getkstar(p1, MassOne, p2, MassTwo);
        if (kstar > ConfOptHighkstarCut)
          continue;

        if (chargeTrack == 1) {
          partSign = 1;
        } else {
          partSign = 1 << 1;
        }

        float invMass;
        if (p2.candidateSelFlag() == 1) {
          invMass = p2.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
        } else {
          invMass = p2.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
        }
        fillFemtoResult(
          invMass,
          p2.pt(),
          p1.pt(),
          p2.bdtBkg(),
          p2.bdtPrompt(),
          p2.bdtFD(),
          kstar,
          FemtoDreamMath::getkT(p1, MassOne, p2, MassTwo),
          FemtoDreamMath::getmT(p1, MassOne, p2, MassTwo),
          collision1.multNtr(),
          collision1.multV0M(),
          partSign,
          processType);

        if (ConfOptUseCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, collision1.magField())) {
            continue;
          }
        }
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        mixedEventCont.setPair<isMc, true>(p1, p2, collision1.multNtr(), collision1.multV0M(), ConfOptUse4D, ConfOptExtendedPlots, ConfOptsmearingByOrigin);
      }
    }
  }

  void processMixedEvent(FilteredColisions const& cols, FilteredFDParticles const& parts, FilteredCharmCands const& /*candidates*/)
  {
    switch (ConfMixingPolicy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent<false>(cols, parts, PartitionTrk1, PartitionHF, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<false>(cols, parts, PartitionTrk1, PartitionHF, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<false>(cols, parts, PartitionTrk1, PartitionHF, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskCharmHadrons, processMixedEvent, "Enable processing mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<femtoDreamPairTaskCharmHadrons>(cfgc)};
}
