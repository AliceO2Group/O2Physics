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

/// \file taskCharmHadronsFemtoDream.cxx.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Ravindra SIngh, GSI, ravindra.singh@cern.ch
/// \author Biao Zhang, Heidelberg University, biao.zhang@cern.ch

#include <vector>
#include <string>

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

struct HfTaskCharmHadronsFemtoDream {

  enum TrackCharge {
    PositiveCharge = 1,
    NegativeCharge = -1
  };

  enum PairSign {
    PairNotDefined = 0,
    LikeSignPair = 1,
    UnLikeSignPair = 2
  };

  /// Binning configurables
  ConfigurableAxis bin4Dkstar{"bin4Dkstar", {1500, 0., 6.}, "binning kstar for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis bin4DMult{"bin4Dmult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis bin4DmT{"bin4DmT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis bin4DmultPercentile{"bin4DmultPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis binInvMass{"binInvMass", {400, 2.10, 2.50}, "InvMass binning"};
  ConfigurableAxis binpTCharm{"binpTCharm", {360, 0, 36}, "pT binning of charm hadron"};
  ConfigurableAxis binTempFitVarTrack{"binTempFitVarTrack", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Track)"};
  ConfigurableAxis binmT{"binmT", {225, 0., 7.5}, "binning mT"};
  ConfigurableAxis binmultTempFit{"binmultTempFit", {1, 0, 1}, "multiplicity Binning for the TempFitVar plot"};
  ConfigurableAxis binMulPercentile{"binMulPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning"};
  ConfigurableAxis binpTTrack{"binpTTrack", {50, 0.5, 10.05}, "pT binning of the pT vs. TempFitVar plot (Track)"};
  ConfigurableAxis binEta{"binEta", {{200, -1.5, 1.5}}, "eta binning"};
  ConfigurableAxis binPhi{"binPhi", {{200, 0, TMath::TwoPi()}}, "phi binning"};
  ConfigurableAxis binkT{"binkT", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis binkstar{"binkstar", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis binNSigmaTPC{"binNSigmaTPC", {1600, -8, 8}, "Binning of Nsigma TPC plot"};
  ConfigurableAxis binNSigmaTOF{"binNSigmaTOF", {3000, -15, 15}, "Binning of the Nsigma TOF plot"};
  ConfigurableAxis binNSigmaTPCTOF{"binNSigmaTPCTOF", {3000, -15, 15}, "Binning of the Nsigma TPC+TOF plot"};
  ConfigurableAxis binTPCClusters{"binTPCClusters", {163, -0.5, 162.5}, "Binning of TPC found clusters plot"};
  Configurable<int> ConfTempFitVarMomentum{"ConfTempFitVarMomentum", 0, "Momentum used for binning: 0 -> pt; 1 -> preco; 2 -> ptpc"};

  /// Particle 2 (Charm Hadrons)
  Configurable<float> charmHadBkgBDTmax{"charmHadBkgBDTmax", 1., "Maximum background bdt score for Charm Hadron (particle 2)"};
  Configurable<int> charmHadCandSel{"charmHadCandSel", 1, "candidate selection for charm hadron"};
  Configurable<int> charmHadMcSel{"charmHadMcSel", 2, "charm hadron selection for mc, partDplusToPiKPi (1), partLcToPKPi (2), partDsToKKPi (4), partXicToPKPi (8)"};
  Configurable<float> charmHadFdBDTmin{"charmHadFdBDTmin", 0., "Minimum feed-down bdt score Charm Hadron (particle 2)"};
  Configurable<float> charmHadFdBDTmax{"charmHadFdBDTmax", 1., "Maximum feed-down bdt score Charm Hadron (particle 2)"};
  Configurable<float> charmHadMaxInvMass{"charmHadMaxInvMass", 2.45, "Maximum invariant mass of Charm Hadron (particle 2)"};
  Configurable<float> charmHadMaxPt{"charmHadMaxPt", 999., "Maximum pT of Charm Hadron (particle 2)"};
  Configurable<float> charmHadMinInvMass{"charmHadMinInvMass", 2.15, "Minimum invariant mass of Charm Hadron (particle 2)"};
  Configurable<float> charmHadMinPt{"charmHadMinPt", 0., "Minimum pT of Charm Hadron (particle 2)"};
  Configurable<int> charmHadPDGCode{"charmHadPDGCode", 4122, "PDG code of particle 2 Charm Hadron"};
  Configurable<float> charmHadPromptBDTmin{"charmHadPromptBDTmin", 0., "Minimum prompt bdt score Charm Hadron (particle 2)"};
  Configurable<float> charmHadPromptBDTmax{"charmHadPromptBDTmax", 1., "Maximum prompt bdt score Charm Hadron (particle 2)"};

  /// General options
  Configurable<float> cprDeltaEtaMax{"cprDeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
  Configurable<float> cprDeltaPhiMax{"cprDeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<bool> cprPlotPerRadii{"cprPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<bool> extendedPlots{"extendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
  Configurable<float> highkstarCut{"highkstarCut", 100000., "Set a cut for high k*, above which the pairs are rejected"};
  Configurable<bool> isMc{"isMc", false, "Set true in the case of a MonteCarlo Run"};
  Configurable<bool> smearingByOrigin{"smearingByOrigin", false, "Obtain the smearing matrix differential in the MC origin of particle 1 and particle 2. High memory consumption. Use with care!"};
  Configurable<bool> use4D{"use4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
  Configurable<bool> useCPR{"useCPR", false, "Close Pair Rejection"};
  ConfigurableAxis dummy{"dummy", {1, 0, 1}, "dummy axis"};

  // Mixing configurables
  ConfigurableAxis mixingBinMult{"mixingBinMult", {VARIABLE_WIDTH, 0.0f, 20.0f, 60.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis mixingBinMultPercentile{"mixingBinMultPercentile", {VARIABLE_WIDTH, 0.0f, 100.f}, "Mixing bins - multiplicity percentile"};
  ConfigurableAxis mixingBinVztx{"mixingBinVztx", {VARIABLE_WIDTH, -10.0f, -4.f, 0.f, 4.f, 10.f}, "Mixing bins - z-vertex"};
  Configurable<int> mixingDepth{"mixingDepth", 5, "Number of events for mixing"};
  Configurable<int> mixingBinPolicy{"mixingBinPolicy", 0, "Binning policy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both"};

  /// Event selection
  struct : ConfigurableGroup {
    std::string prefix = "eventSel";
    Configurable<int> multMin{"multMin", 0, "Minimum Multiplicity (MultNtr)"};
    Configurable<int> multMax{"multMax", 99999, "Maximum Multiplicity (MultNtr)"};
    Configurable<float> multPercentileMin{"multPercentileMin", 0, "Maximum Multiplicity Percentile"};
    Configurable<float> multPercentileMax{"multPercentileMax", 100, "Minimum Multiplicity Percentile"};
  } eventSel;

  /// Particle 1 (track)
  Configurable<femtodreamparticle::cutContainerType> cutBitTrack1{"cutBitTrack1", 5542474, "Particle 1 (Track) - Selection bit from cutCulator"};
  Configurable<int> pdgCodeTrack1{"pdgCodeTrack1", 2212, "PDG code of Particle 1 (Track)"};
  Configurable<float> pidThresTrack1{"pidThresTrack1", 0.75, "Momentum threshold for PID selection for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> tpcBitTrack1{"tpcBitTrack1", 4, "PID TPC bit from cutCulator for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> tpcTofBitTrack1{"tpcTofBitTrack1", 2, "PID TPCTOF bit from cutCulator for particle 1 (Track)"};
  Configurable<float> etaTrack1Max{"etaTrack1Max", 10., "Maximum eta of partricle 1 (Track)"};
  Configurable<float> ptTrack1Max{"ptTrack1Max", 999., "Maximum pT of partricle 1 (Track)"};
  Configurable<float> etaTrack1Min{"etaTrack1Min", -10., "Minimum eta of partricle 1 (Track)"};
  Configurable<float> ptTrack1Min{"ptTrack1Min", 0., "Minimum pT of partricle 1 (Track)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{mixingBinVztx, mixingBinMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{mixingBinVztx, mixingBinMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{mixingBinVztx, mixingBinMult, mixingBinMultPercentile}, true};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron> pairCloseRejectionSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron> pairCloseRejectionME;
  Filter eventMultiplicity = aod::femtodreamcollision::multNtr >= eventSel.multMin && aod::femtodreamcollision::multNtr <= eventSel.multMax;
  Filter eventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= eventSel.multPercentileMin && aod::femtodreamcollision::multV0M <= eventSel.multPercentileMax;
  Filter hfCandSelFilter = aod::fdhf::candidateSelFlag >= static_cast<int8_t>(charmHadCandSel.value);
  Filter hfMcSelFilter = nabs(aod::fdhf::flagMc) == static_cast<int8_t>(charmHadMcSel.value);
  Filter trackEtaFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta < etaTrack1Max, true);
  Filter trackEtaFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta > etaTrack1Min, true);
  Filter trackPtFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt < ptTrack1Max, true);
  Filter trackPtFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt > ptTrack1Min, true);

  using FilteredCharmCands = soa::Filtered<aod::FDHfCand>;
  using FilteredCharmCand = FilteredCharmCands::iterator;

  using FilteredCharmMcCands = soa::Filtered<soa::Join<aod::FDHfCand, aod::FDHfCandMC>>;
  using FilteredCharmMcCand = FilteredCharmMcCands::iterator;

  using FilteredColisions = soa::Filtered<soa::Join<FDCollisions, FDColMasks>>;
  using FilteredColision = FilteredColisions::iterator;

  using FilteredMcColisions = soa::Filtered<soa::Join<aod::FDCollisions, FDColMasks, aod::FDMCCollLabels>>;
  using FilteredMcColision = FilteredMcColisions::iterator;

  using FilteredFDMcParts = soa::Filtered<soa::Join<aod::FDParticles, aod::FDParticlesIndex, aod::FDExtParticles, aod::FDMCLabels, aod::FDExtMCLabels>>;
  using FilteredFDMcPart = FilteredFDMcParts::iterator;

  using FilteredFDParticles = soa::Filtered<soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDParticlesIndex>>;
  using FilteredFDParticle = FilteredFDParticles::iterator;

  femtodreamcollision::BitMaskType BitMask = 1 << 0;

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> allTrackHisto;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 5> selectedTrackHisto;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;
  /// Histogram output
  HistogramRegistry registry{"CorrelationsAndQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryMixQa{"registryMixQa"};
  HistogramRegistry registryCharmHadronQa{"registryCharmHadronQa"};
  /// Partition for particle 1

  Partition<FilteredFDParticles> partitionTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && (ncheckbit(aod::femtodreamparticle::cut, cutBitTrack1)) && ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= pidThresTrack1, ncheckbit(aod::femtodreamparticle::pidcut, tpcBitTrack1), ncheckbit(aod::femtodreamparticle::pidcut, tpcTofBitTrack1));

  Partition<FilteredFDMcParts> partitionMcTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 (ncheckbit(aod::femtodreamparticle::cut, cutBitTrack1)) &&
                                                 ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= pidThresTrack1, ncheckbit(aod::femtodreamparticle::pidcut, tpcBitTrack1), ncheckbit(aod::femtodreamparticle::pidcut, tpcTofBitTrack1));

  /// Partition for particle 2
  Partition<FilteredCharmCands> partitionCharmHadron = aod::fdhf::bdtBkg < charmHadBkgBDTmax && aod::fdhf::bdtFD < charmHadFdBDTmax && aod::fdhf::bdtFD > charmHadFdBDTmin&& aod::fdhf::bdtPrompt<charmHadPromptBDTmax && aod::fdhf::bdtPrompt> charmHadPromptBDTmin;
  Partition<FilteredCharmMcCands> partitionMcCharmHadron = aod::fdhf::originMcRec == 1 || aod::fdhf::originMcRec == 2;

  float massOne = o2::analysis::femtoDream::getMass(pdgCodeTrack1);
  float massTwo = o2::analysis::femtoDream::getMass(charmHadPDGCode);
  int8_t partSign = 0;
  int64_t processType = 0;

  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  Produces<o2::aod::FDResultsHF> fillFemtoResult;

  void init(InitContext& /*context*/)
  {
    // setup columnpolicy for binning
    colBinningMult = {{mixingBinVztx, mixingBinMult}, true};
    colBinningMultPercentile = {{mixingBinVztx, mixingBinMultPercentile}, true};
    colBinningMultMultPercentile = {{mixingBinVztx, mixingBinMult, mixingBinMultPercentile}, true};
    eventHisto.init(&registry);
    allTrackHisto.init(&registry, binmultTempFit, binMulPercentile, binpTTrack, binEta, binPhi, binTempFitVarTrack, binNSigmaTPC, binNSigmaTOF, binNSigmaTPCTOF, binTPCClusters, dummy, isMc, pdgCodeTrack1, true);
    selectedTrackHisto.init(&registry, binmultTempFit, binMulPercentile, binpTTrack, binEta, binPhi, binTempFitVarTrack, binNSigmaTPC, binNSigmaTOF, binNSigmaTPCTOF, binTPCClusters, dummy, isMc, pdgCodeTrack1, true);

    sameEventCont.init<true>(&registry,
                             binkstar, binpTTrack, binkT, binmT, mixingBinMult, mixingBinMultPercentile,
                             bin4Dkstar, bin4DmT, bin4DMult, bin4DmultPercentile,
                             isMc, use4D, extendedPlots,
                             highkstarCut,
                             smearingByOrigin, binInvMass);

    sameEventCont.setPDGCodes(pdgCodeTrack1, charmHadPDGCode);
    mixedEventCont.init<true>(&registry,
                              binkstar, binpTTrack, binkT, binmT, mixingBinMult, mixingBinMultPercentile,
                              bin4Dkstar, bin4DmT, bin4DMult, bin4DmultPercentile,
                              isMc, use4D, extendedPlots,
                              highkstarCut,
                              smearingByOrigin, binInvMass);

    mixedEventCont.setPDGCodes(pdgCodeTrack1, charmHadPDGCode);
    registryMixQa.add("MixingQA/hSECollisionBins", "; bin; Entries", kTH1F, {{120, -0.5, 119.5}});
    registryMixQa.add("MixingQA/hSECollisionPool", "; Vz (cm); Mul", kTH2F, {{100, -10, 10}, {200, 0, 200}});
    registryMixQa.add("MixingQA/hMECollisionBins", "; bin; Entries", kTH1F, {{120, -0.5, 119.5}});
    registryCharmHadronQa.add("CharmHadronQA/hPtVsMass", "; #it{p}_{T} (GeV/#it{c}); inv. mass (GeV/#it{c}^{2})", kTH2F, {binpTCharm, binInvMass});

    pairCleaner.init(&registry);
    if (useCPR.value) {
      pairCloseRejectionSE.init(&registry, &registry, cprDeltaPhiMax.value, cprDeltaEtaMax.value, cprPlotPerRadii.value, 1);
      pairCloseRejectionME.init(&registry, &registry, cprDeltaPhiMax.value, cprDeltaEtaMax.value, cprPlotPerRadii.value, 2);
    }
  }

  template <typename CollisionType>
  void fillCollision(CollisionType const& col)
  {
    registryMixQa.fill(HIST("MixingQA/hSECollisionBins"), colBinningMult.getBin({col.posZ(), col.multNtr()}));
    registryMixQa.fill(HIST("MixingQA/hSECollisionPool"), col.posZ(), col.multNtr());
  }

  /// This function processes the same event and takes care of all the histogramming
  template <bool isMc, typename PartitionType, typename CandType, typename TableTracks, typename Collision>
  void doSameEvent(PartitionType& sliceTrk1, CandType& sliceCharmHad, TableTracks const& parts, Collision const& col)
  {
    fillCollision(col);

    processType = 1; // for same event

    for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceTrk1, sliceCharmHad))) {

      if (p1.trackId() == p2.prong0Id() || p1.trackId() == p2.prong1Id() || p1.trackId() == p2.prong2Id())
        continue;

      if (useCPR.value) {
        if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
          continue;
        }
      }

      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }

      // proton track charge
      float chargeTrack = 0.;
      if ((p1.cut() & 2) == 2) {
        chargeTrack = PositiveCharge;
      } else {
        chargeTrack = NegativeCharge;
      }

      int pairSign = 0;
      if (chargeTrack == p2.charge()) {
        pairSign = LikeSignPair;
      } else {
        pairSign = UnLikeSignPair;
      }

      float kstar = FemtoDreamMath::getkstar(p1, massOne, p2, massTwo);
      if (kstar > highkstarCut) {
        continue;
      }

      // if (chargeTrack == 1) {
      //   partSign = 1;
      // } else {
      //   partSign = 1 << 1;
      // }

      float invMass;
      if (p2.candidateSelFlag() == 1) {
        invMass = p2.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
      } else {
        invMass = p2.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
      }

      if (invMass < charmHadMinInvMass || invMass > charmHadMaxInvMass) {
        continue;
      }

      if (p2.pt() < charmHadMinPt || p2.pt() > charmHadMaxPt) {
        continue;
      }
      /// Filling QA histograms of the selected tracks
      selectedTrackHisto.fillQA<isMc, true>(p1, static_cast<aod::femtodreamparticle::MomentumType>(ConfTempFitVarMomentum.value), col.multNtr(), col.multV0M());

      int charmHadMc = 0;
      int originType = 0;
      if constexpr (isMc) {
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
        FemtoDreamMath::getkT(p1, massOne, p2, massTwo),
        FemtoDreamMath::getmT(p1, massOne, p2, massTwo),
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

    // Mixed events that contain the pair of interest
    processType = 2; // for mixed event

    Partition<CollisionType> PartitionMaskedCol1 = (aod::femtodreamcollision::bitmaskTrackOne & BitMask) == BitMask;
    PartitionMaskedCol1.bindTable(cols);

    Partition<CollisionType> PartitionMaskedCol2 = (aod::femtodreamcollision::bitmaskTrackTwo & BitMask) == BitMask;
    PartitionMaskedCol2.bindTable(cols);

    for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockFullIndexPolicy(policy, mixingDepth.value, -1, *PartitionMaskedCol1.mFiltered, *PartitionMaskedCol2.mFiltered))) {
      // make sure that tracks in the same events are not mixed
      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      const int multiplicityCol = collision1.multNtr();
      registryMixQa.fill(HIST("MixingQA/hMECollisionBins"), colBinningMult.getBin({collision1.posZ(), multiplicityCol}));

      auto sliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto sliceCharmHad = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceTrk1, sliceCharmHad))) {

        if (useCPR.value) {
          if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
            continue;
          }
        }
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }

        float chargeTrack = 0.;
        if ((p1.cut() & 2) == 2) {
          chargeTrack = PositiveCharge;
        } else {
          chargeTrack = NegativeCharge;
        }

        int pairSign = 0;
        if (chargeTrack == p2.charge()) {
          pairSign = LikeSignPair;
        } else {
          pairSign = UnLikeSignPair;
        }

        float kstar = FemtoDreamMath::getkstar(p1, massOne, p2, massTwo);
        if (kstar > highkstarCut) {
          continue;
        }
        float invMass;
        if (p2.candidateSelFlag() == 1) {
          invMass = p2.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
        } else {
          invMass = p2.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
        }

        if (invMass < charmHadMinInvMass || invMass > charmHadMaxInvMass) {
          continue;
        }

        if (p2.pt() < charmHadMinPt || p2.pt() > charmHadMaxPt) {
          continue;
        }

        int charmHadMc = 0;
        int originType = 0;
        if constexpr (isMc) {
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
          FemtoDreamMath::getkT(p1, massOne, p2, massTwo),
          FemtoDreamMath::getmT(p1, massOne, p2, massTwo),
          collision1.multNtr(),
          collision1.multV0M(),
          p2.charge(),
          pairSign,
          processType,
          charmHadMc,
          originType);

        // if constexpr (!isMc) mixedEventCont.setPair<isMc, true>(p1, p2, collision1.multNtr(), collision1.multV0M(), use4D, extendedPlots, smearingByOrigin);
        mixedEventCont.setPair<isMc, true>(p1, p2, collision1.multNtr(), collision1.multV0M(), use4D, extendedPlots, smearingByOrigin);
      }
    }
  }

  void processSameEvent(FilteredColision const& col,
                        FilteredFDParticles const& parts,
                        FilteredCharmCands const&)
  {
    eventHisto.fillQA(col);
    auto sliceTrk1 = partitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto sliceCharmHad = partitionCharmHadron->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    /// Filling QA histograms of the all tracks and all charm hadrons before pairing
    for (auto const& part : sliceTrk1) {
      allTrackHisto.fillQA<false, true>(part, static_cast<aod::femtodreamparticle::MomentumType>(ConfTempFitVarMomentum.value), col.multNtr(), col.multV0M());
    }
    for (auto const& part : sliceCharmHad) {
      float invMass;
      if (part.candidateSelFlag() == 1) {
        invMass = part.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
      } else {
        invMass = part.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
      }
      registryCharmHadronQa.fill(HIST("CharmHadronQA/hPtVsMass"), part.pt(), invMass);
    }

    if ((col.bitmaskTrackOne() & BitMask) != BitMask || (col.bitmaskTrackTwo() & BitMask) != BitMask) {
      return;
    }
    doSameEvent<false>(sliceTrk1, sliceCharmHad, parts, col);
  }
  PROCESS_SWITCH(HfTaskCharmHadronsFemtoDream, processSameEvent, "Enable processing same event", false);

  void processMixedEvent(FilteredColisions const& cols,
                         FilteredFDParticles const& parts,
                         FilteredCharmCands const&)
  {
    switch (mixingBinPolicy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent<false>(cols, parts, partitionTrk1, partitionCharmHadron, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<false>(cols, parts, partitionTrk1, partitionCharmHadron, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<false>(cols, parts, partitionTrk1, partitionCharmHadron, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsFemtoDream, processMixedEvent, "Enable processing mixed events", false);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMc(FilteredMcColision const& col,
                          FilteredFDMcParts const& parts,
                          o2::aod::FDMCParticles const&,
                          o2::aod::FDExtMCParticles const&,
                          FilteredCharmMcCands const&)
  {
    eventHisto.fillQA(col);

    auto sliceMcTrk1 = partitionMcTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto sliceMcCharmHad = partitionMcCharmHadron->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    if (sliceMcTrk1.size() == 0 && sliceMcCharmHad.size() == 0) {
      return;
    }
    /// Filling QA histograms of the all mc tracks before pairing
    for (auto const& part : sliceMcTrk1) {
      allTrackHisto.fillQA<true, true>(part, static_cast<aod::femtodreamparticle::MomentumType>(ConfTempFitVarMomentum.value), col.multNtr(), col.multV0M());
    }

    if ((col.bitmaskTrackOne() & BitMask) != BitMask || (col.bitmaskTrackTwo() & BitMask) != BitMask) {
      return;
    }
    doSameEvent<true>(sliceMcTrk1, sliceMcCharmHad, parts, col);
  }
  PROCESS_SWITCH(HfTaskCharmHadronsFemtoDream, processSameEventMc, "Enable processing same event for Monte Carlo", false);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMc(FilteredMcColisions const& cols,
                           FilteredFDMcParts const& parts,
                           o2::aod::FDMCParticles const&,
                           o2::aod::FDExtMCParticles const&,
                           FilteredCharmMcCands const&)
  {
    switch (mixingBinPolicy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent<true>(cols, parts, partitionMcTrk1, partitionMcCharmHadron, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<true>(cols, parts, partitionMcTrk1, partitionMcCharmHadron, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<true>(cols, parts, partitionMcTrk1, partitionMcCharmHadron, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsFemtoDream, processMixedEventMc, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCharmHadronsFemtoDream>(cfgc)};
}
