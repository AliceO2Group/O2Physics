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

  enum TrackCharge {
    positiveCharge = 1,
    negativeCharge = -1
  };

  enum PairSign {
    pairNotDefined = 0,
    likeSignPair = 1,
    unLikeSignPair = 2
  };

  /// General options
  Configurable<bool> CPRPlotPerRadii{"CPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<bool> useCPR{"useCPR", false, "Close Pair Rejection"};
  Configurable<float> CPRdeltaEtaMax{"CPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
  Configurable<float> CPRdeltaPhiMax{"CPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<bool> extendedPlots{"extendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
  Configurable<float> highkstarCut{"highkstarCut", 100000., "Set a cut for high k*, above which the pairs are rejected"};
  //Configurable<bool> isMc{"isMc", true, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<bool> isMc{"isMc", false, "Set true in the case of a MonteCarlo Run"};
  Configurable<bool> smearingByOrigin{"smearingByOrigin", false, "Obtain the smearing matrix differential in the MC origin of particle 1 and particle 2. High memory consumption. Use with care!"};
  Configurable<bool> use4D{"use4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
  ConfigurableAxis dummy{"dummy", {1, 0, 1}, "dummy axis"};

  /// Event selection
  //Configurable<int> event_maxMult{"event_maxMult", 99999, "Maximum Multiplicity (MultNtr)"};
  //Configurable<float> event_maxMultPercentile{"event_maxMultPercentile", 100, "Maximum Multiplicity Percentile"};
  //Configurable<int> event_minMult{"event_minMult", 0, "Minimum Multiplicity (MultNtr)"};
  //Configurable<float> event_minMultPercentile{"event_minMultPercentile", 0, "Minimum Multiplicity Percentile"};

  /// Particle 1 (track)
  Configurable<femtodreamparticle::cutContainerType> trk1_CutBit{"trk1_CutBit", 5542474, "Particle 1 (Track) - Selection bit from cutCulator"};
  Configurable<int> trk1_PDGCode{"trk1_PDGCode", 2212, "PDG code of Particle 1 (Track)"};
  Configurable<float> trk1_PIDThres{"trk1_PIDThres", 0.75, "Momentum threshold for PID selection for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> trk1_TPCBit{"trk1_TPCBit", 4, "PID TPC bit from cutCulator for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> trk1_TPCTOFBit{"trk1_TPCTOFBit", 2, "PID TPCTOF bit from cutCulator for particle 1 (Track)"};
  Configurable<float> trk1_maxEta{"trk1_maxEta", 10., "Maximum eta of partricle 1 (Track)"};
  Configurable<float> trk1_maxPt{"trk1_maxPt", 999., "Maximum pT of partricle 1 (Track)"};
  Configurable<float> trk1_minEta{"trk1_minEta", -10., "Minimum eta of partricle 1 (Track)"};
  Configurable<float> trk1_minPt{"trk1_minPt", 0., "Minimum pT of partricle 1 (Track)"};

  /// Particle 2 (Charm Hadrons)
  Configurable<float> charmHad_bkgBDT{"charmHad_bkgBDT", 1., "Maximum background bdt score for Charm Hadron (particle 2)"};
  Configurable<int8_t> charmHad_CandSel{"charmHad_CandSel", 1, "candidate selection for charm hadron"};
  Configurable<int8_t> charmHad_McSel{"charmHad_McSel", 2, "charm hadron selection for mc, partDplusToPiKPi (1), partLcToPKPi (2), partDsToKKPi (4), partXicToPKPi (8)"};
  Configurable<float> charmHad_fdBDT{"charmHad_fdBDT", 0., "Minimum feed-down bdt score Charm Hadron (particle 2)"};
  Configurable<float> charmHad_maxInvMass{"charmHad_maxInvMass", 2.45, "Maximum invariant mass of Charm Hadron (particle 2)"};
  Configurable<float> charmHad_maxPt{"charmHad_maxPt", 999., "Maximum pT of Charm Hadron (particle 2)"};
  Configurable<float> charmHad_minInvMass{"charmHad_minInvMass", 2.15, "Minimum invariant mass of Charm Hadron (particle 2)"};
  Configurable<float> charmHad_minPt{"charmHad_minPt", 0., "Minimum pT of Charm Hadron (particle 2)"};
  Configurable<int> charmHad_PDGCode{"charmHad_PDGCode", 4122, "PDG code of particle 2 Charm Hadron"};
  Configurable<float> charmHad_promptBDT{"charmHad_promptBDT", 0., "Minimum prompt bdt score Charm Hadron (particle 2)"};

  /// Binning configurables
  ConfigurableAxis bin4Dkstar{"bin4Dkstar", {1500, 0., 6.}, "binning kstar for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis bin4DMult{"bin4Dmult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis bin4DmT{"bin4DmT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis bin4DmultPercentile{"bin4DmultPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis binInvMass{"binInvMass", {300, 2.15, 2.45}, "InvMass binning"};
  //ConfigurableAxis binTempFitVarHF{"binTempFitVarHF", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0)"};
  //ConfigurableAxis binTempFitVarHFChild{"binTempFitVarHFChild", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0 child)"};
  ConfigurableAxis binTempFitVarTrack{"binTempFitVarTrack", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Track)"};
  ConfigurableAxis binmT{"binmT", {225, 0., 7.5}, "binning mT"};
  ConfigurableAxis binmultTempFit{"binmultTempFit", {1, 0, 1}, "multiplicity Binning for the TempFitVar plot"};
  ConfigurableAxis binpT{"binpT", {20, 0.5, 4.05}, "pT binning"};
  //ConfigurableAxis binpTHF{"binpTHF", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0)"};
  //ConfigurableAxis binpTHFChild{"binpTHFChild", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis binpTTrack{"binpTTrack", {50, 0.5, 10.05}, "pT binning of the pT vs. TempFitVar plot (Track)"};
  ConfigurableAxis binkT{"binkT", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis binkstar{"binkstar", {1500, 0., 6.}, "binning kstar"};

  // Mixing configurables
  ConfigurableAxis mixingBinMult{"mixingBinMult", {VARIABLE_WIDTH, 0.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis mixingBinMultPercentile{"mixingBinMultPercentile", {VARIABLE_WIDTH, 0.0f, 100.f}, "Mixing bins - multiplicity percentile"};
  ConfigurableAxis mixingBinVztx{"mixingBinVztx", {VARIABLE_WIDTH, -10.0f, -4.f, 0.f, 4.f, 10.f}, "Mixing bins - z-vertex"};
  Configurable<int> mixingDepth{"mixingDepth", 5, "Number of events for mixing"};
  Configurable<int> mixingPolicy{"mixingBinPolicy", 0, "Binning policy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both"};

  /// Event selection
  struct : ConfigurableGroup {
    Configurable<int> MultMin{"EventSel.MultMin", 0, "Minimum Multiplicity (MultNtr)"};
    Configurable<int> MultMax{"EventSel.MultMax", 99999, "Maximum Multiplicity (MultNtr)"};
    Configurable<float> MultPercentileMin{"EventSel.MultPercentileMin", 0, "Maximum Multiplicity Percentile"};
    Configurable<float> MultPercentileMax{"EventSel.MultPercentileMax", 100, "Minimum Multiplicity Percentile"};
  } EventSel;

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{mixingBinVztx, mixingBinMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{mixingBinVztx, mixingBinMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{mixingBinVztx, mixingBinMult, mixingBinMultPercentile}, true};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron> pairCloseRejection;
  
  Filter EventMultiplicity = aod::femtodreamcollision::multNtr >= EventSel.MultMin && aod::femtodreamcollision::multNtr <= EventSel.MultMax;
  Filter EventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= EventSel.MultPercentileMin && aod::femtodreamcollision::multV0M <= EventSel.MultPercentileMax;
  Filter trackEtaFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta < trk1_maxEta, true);
  Filter trackEtaFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta > trk1_minEta, true);
  Filter trackPtFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt < trk1_maxPt, true);
  Filter trackPtFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt > trk1_minPt, true);
  Filter hfCandSelFilter = aod::fdhf::candidateSelFlag >= charmHad_CandSel.value;
  Filter hfMcSelFilter = nabs(aod::fdhf::flagMc) == charmHad_McSel.value;

  using FilteredCharmCands = soa::Filtered<aod::FDHfCand>;
  using FilteredCharmCand = FilteredCharmCands::iterator;

  using FilteredCharmMcCands = soa::Filtered<soa::Join<aod::FDHfCand, aod::FDHfCandMC>>;
  using FilteredCharmMcCand = FilteredCharmMcCands::iterator;

  using FilteredColisions = FDCollisions;
  using FilteredColision = FilteredColisions::iterator;

  using FilteredMcColisions = soa::Filtered<soa::Join<aod::FDCollisions, aod::FDMCCollLabels>>;
  using FilteredMcColision = FilteredMcColisions::iterator;

  using FilteredFDMcParts = soa::Filtered<soa::Join<aod::FDParticles, aod::FDParticlesIndex, aod::FDMCLabels>>;
  using FilteredFDMcPart = FilteredFDMcParts::iterator;

  using FilteredFDParticles = soa::Filtered<soa::Join<aod::FDParticles, aod::FDParticlesIndex>>;
  using FilteredFDParticle = FilteredFDParticles::iterator;

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Partition for particle 1

  Partition<FilteredFDParticles> PartitionTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack));

  Partition<FilteredFDMcParts> PartitionMcTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 (ncheckbit(aod::femtodreamparticle::cut, trk1_CutBit)) &&
                                                 ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= trk1_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, trk1_TPCBit), ncheckbit(aod::femtodreamparticle::pidcut, trk1_TPCTOFBit));

  /// Partition for particle 2
  Partition<FilteredCharmCands> PartitionHF = aod::fdhf::bdtBkg < charmHad_bkgBDT && aod::fdhf::bdtPrompt > charmHad_promptBDT;
  Partition<FilteredCharmMcCands> PartitionMcHF = aod::fdhf::originMcRec == 0 || aod::fdhf::originMcRec == 1 || aod::fdhf::originMcRec == 2;
  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"CorrelationsHF", {}, OutputObjHandlingPolicy::AnalysisObject};

  float MassOne = o2::analysis::femtoDream::getMass(trk1_PDGCode);
  float MassTwo = o2::analysis::femtoDream::getMass(charmHad_PDGCode);
  int8_t partSign = 0;
  int64_t processType = 0;

  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  Produces<o2::aod::FDResultsHF> fillFemtoResult;

  void init(InitContext& context)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry, binmultTempFit, dummy, binpTTrack, dummy, dummy, binTempFitVarTrack, dummy, dummy, dummy, dummy, dummy, isMc, trk1_PDGCode);

    sameEventCont.init<true>(&resultRegistry,
                             binkstar, binpT, binkT, binmT, mixingBinMult, mixingBinMultPercentile,
                             bin4Dkstar, bin4DmT, bin4DMult, bin4DmultPercentile,
                             isMc, use4D, extendedPlots,
                             highkstarCut,
                             smearingByOrigin, binInvMass);

    sameEventCont.setPDGCodes(trk1_PDGCode, charmHad_PDGCode);
    mixedEventCont.init<true>(&resultRegistry,
                              binkstar, binpT, binkT, binmT, mixingBinMult, mixingBinMultPercentile,
                              bin4Dkstar, bin4DmT, bin4DMult, bin4DmultPercentile,
                              isMc, use4D, extendedPlots,
                              highkstarCut,
                              smearingByOrigin, binInvMass);

    mixedEventCont.setPDGCodes(trk1_PDGCode, charmHad_PDGCode);
    pairCleaner.init(&qaRegistry);
    if (useCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, CPRdeltaPhiMax.value, CPRdeltaEtaMax.value, CPRPlotPerRadii.value);
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  template <bool isMc, typename PartitionType, typename CandType, typename TableTracks, typename TableCandidates, typename Collision>
  void doSameEvent(PartitionType& SliceTrk1, CandType& SliceCharmHad, TableTracks const& parts, TableCandidates const& candidates, Collision const& col)
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

      int pairSign = 0;
      if (chargeTrack == p2.charge()) {
        pairSign = likeSignPair;
      } else { 
        pairSign = unLikeSignPair;
      }

      float kstar = FemtoDreamMath::getkstar(p1, MassOne, p2, MassTwo);
      if (kstar > highkstarCut)
        continue;

     // if (chargeTrack == 1) {
     //   partSign = 1;
     // } else {
     //   partSign = 1 << 1;
     // }

      if (useCPR.value) {
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

      if (invMass < charmHad_minInvMass || invMass > charmHad_maxInvMass)
        continue;
      if (p2.pt() < charmHad_minPt || p2.pt() > charmHad_maxPt)
        continue;
      int charmHadMc = 0;
      int originType = 0;
      if constexpr (isMc){
      charmHadMc = p2.flagMc();
      originType = p2.originMcRec();
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
        col.multNtr(),
        col.multV0M(),
        p2.charge(),
        pairSign,
        processType,
        charmHadMc,
        originType);

      sameEventCont.setPair<isMc, true>(p1, p2, col.multNtr(), col.multV0M(), use4D, extendedPlots, smearingByOrigin);
    }
  }

  template <bool isMc, typename CollisionType, typename PartType, typename PartitionType1, typename PartitionType2, typename BinningType>
  void doMixedEvent(CollisionType const& cols, PartType const& parts, PartitionType1& part1, PartitionType2& part2, BinningType policy)
  {
    processType = 1 << 1; // for mixed event

    for (auto const& [collision1, collision2] : soa::selfCombinations(policy, mixingDepth.value, -1, cols, cols)) {

      auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceCharmHad = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceCharmHad))) {

        float chargeTrack = 0.;
      if ((p1.cut() & 1) == 1) {
        chargeTrack = positiveCharge;
      } else {
        chargeTrack = negativeCharge;
      }

      int pairSign = 0;
      if (chargeTrack == p2.charge()) {
        pairSign = likeSignPair;
      } else { 
        pairSign = unLikeSignPair;
      }

      float kstar = FemtoDreamMath::getkstar(p1, MassOne, p2, MassTwo);
      if (kstar > highkstarCut) continue;

        float invMass;
        if (p2.candidateSelFlag() == 1) {
          invMass = p2.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
        } else {
          invMass = p2.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
        }

      int charmHadMc = 0;
      int originType = 0;
      if constexpr (isMc){
      charmHadMc = p2.flagMc();
      originType = p2.originMcRec();
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
        p2.charge(),
        pairSign,
        processType,
        charmHadMc,
        originType);

        if (useCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, collision1.magField())) {
            continue;
          }
        }
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        //if constexpr (!isMc) mixedEventCont.setPair<isMc, true>(p1, p2, collision1.multNtr(), collision1.multV0M(), use4D, extendedPlots, smearingByOrigin);
        mixedEventCont.setPair<isMc, true>(p1, p2, collision1.multNtr(), collision1.multV0M(), use4D, extendedPlots, smearingByOrigin);
      }
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

  void processMixedEvent(FilteredColisions const& cols, FilteredFDParticles const& parts, FilteredCharmCands const& candidates)
  {
    switch (mixingPolicy.value) {
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

/// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMc(FilteredMcColision& col,
                          o2::aod::FDMCCollisions&,
                          FilteredFDMcParts const& parts,
                          o2::aod::FDMCParticles const&, FilteredCharmMcCands const& candidates)
  {
    auto SliceMcTrk1 = PartitionMcTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceMcCharmHad = PartitionMcHF->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    
    if (SliceMcTrk1.size() == 0 && SliceMcCharmHad.size() == 0) {
      return;
    }
    doSameEvent<true>(SliceMcTrk1, SliceMcCharmHad, parts, candidates, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskCharmHadrons, processSameEventMc, "Enable processing same event for Monte Carlo", false);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMc(FilteredMcColisions const& cols, o2::aod::FDMCCollisions const&, FilteredFDMcParts const& parts, o2::aod::FDMCParticles const&, FilteredCharmMcCands const& candidates)
  {
    switch (mixingPolicy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent<true>(cols, parts, PartitionMcTrk1, PartitionMcHF, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<true>(cols, parts, PartitionMcTrk1, PartitionMcHF, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<true>(cols, parts, PartitionMcTrk1, PartitionMcHF, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskCharmHadrons, processMixedEventMc, "Enable processing mixed events MC", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<femtoDreamPairTaskCharmHadrons>(cfgc)};
}
