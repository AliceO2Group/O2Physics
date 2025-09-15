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

/// \file femtoDreamPairCascadeCascade.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two cascades
/// \author Andi Mathis, Anton Riedel, Georgios Mantzaridis, Oton Vazquez Doce.

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamContainer.h"
#include "PWGCF/FemtoDream/Core/femtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include <sys/stat.h>

#include <cstdint>
#include <string>
#include <vector>
using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;
struct FemtoDreamPairCascadeCascade {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  /// General options
  struct : ConfigurableGroup {
    std::string prefix = std::string("Option");
    Configurable<bool> sameSpecies{"sameSpecies", false, "Set to true if particle 1 and particle 2 are the same species"};
    Configurable<bool> isMC{"isMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<bool> use4D{"use4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
    Configurable<bool> extendedPlots{"extendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
    Configurable<float> highkstarCut{"highkstarCut", -1., "Set a cut for high k*, above which the pairs are rejected. Set it to -1 to deactivate it"};
    Configurable<bool> cprOn{"cprOn", true, "Close Pair Rejection"};
    Configurable<bool> cprOld{"cprOld", false, "Set to FALSE to use fixed version of CPR (for testing now, will be default soon)"};
    Configurable<bool> cprPlotPerRadii{"cprPlotPerRadii", false, "Plot CPR per radii"};
    Configurable<float> cprDeltaPhiMax{"cprDeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
    Configurable<float> cprDeltaEtaMax{"cprDeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
    Configurable<bool> dcaCutPtDep{"dcaCutPtDep", false, "Use pt dependent dca cut"};
    Configurable<bool> mixEventWithPairs{"mixEventWithPairs", false, "Only use events that contain particle 1 and partile 2 for the event mixing"};
    Configurable<bool> smearingByOrigin{"smearingByOrigin", false, "Obtain the smearing matrix differential in the MC origin of particle 1 and particle 2. High memory consumption. Use with care!"};
    ConfigurableAxis dummy{"dummy", {1, 0, 1}, "dummy axis"};
  } Option;

  /// Event selection
  struct : ConfigurableGroup {
    std::string prefix = std::string("EventSel");
    Configurable<int> multMin{"multMin", 0, "Minimum Multiplicity (MultNtr)"};
    Configurable<int> multMax{"multMax", 99999, "Maximum Multiplicity (MultNtr)"};
    Configurable<float> multPercentileMin{"multPercentileMin", 0, "Minimum Multiplicity Percentile"};
    Configurable<float> multPercentileMax{"multPercentileMax", 100, "Maximum Multiplicity Percentile"};
  } EventSel;

  // Filter EventMultiplicity = aod::femtodreamcollision::multNtr >= EventSel.multMin && aod::femtodreamcollision::multNtr <= EventSel.multMax;
  // Filter EventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= EventSel.multPercentileMin && aod::femtodreamcollision::multV0M <= EventSel.multPercentileMax;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;
  // using FilteredCollisions = soa::Filtered<FDCollisions>;
  using FilteredCollisions = FDCollisions;
  using FilteredCollision = FilteredCollisions::iterator;
  using FDMCParts = soa::Join<aod::FDParticles, aod::FDMCLabels>;
  using FDMCPart = FDMCParts::iterator;
  femtodreamcollision::BitMaskType bitMask = 1; //???????????????????

  /// Cascade 1 (Cascade)
  struct : ConfigurableGroup {
    std::string prefix = std::string("Cascade1");
    Configurable<int> pdgCode{"pdgCode", 3334, "PDG code of Particle 1 (Cascade)"};
    Configurable<femtodreamparticle::cutContainerType> cutBit{"cutBit", 5542474, "Particle 1 (Cascade) - Selection bit from cutCulator"};
    Configurable<femtodreamparticle::cutContainerType> childPosCutBit{"childPosCutBit", 278, "Selection bit for positive child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> childPosTPCBit{"childPosTPCBit", 1024, "PID TPC bit for positive child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> childNegCutBit{"childNegCutBit", 277, "Selection bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> childNegTPCBit{"childNegTPCBit", 4096, "PID TPC bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> childBachCutBit{"childBachCutBit", 277, "Selection bit for bachelor child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> childBachTPCBit{"childBachTPCBit", 64, "PID TPC bit for bachelor child of Cascade"};

    Configurable<float> invMassMin{"invMassMin", 1.6, "Minimum invariant mass of Partricle 1 (Cascade)"};
    Configurable<float> invMassMax{"invMassMax", 1.8, "Maximum invariant mass of Partricle 1 (Cascade)"};
    Configurable<float> invMassV0DaughMin{"invMassV0DaughMin", 0., "Minimum invariant mass of the V0 Daughter"};
    Configurable<float> invMassV0DaughMax{"invMassV0DaughMax", 999., "Maximum invariant mass of the V0 Daughter"};
    Configurable<float> ptMin{"ptMin", 0., "Minimum pT of Particle 2 (Cascade)"};
    Configurable<float> ptMax{"ptMax", 999., "Maximum pT of Particle 2 (Cascade)"};
    Configurable<float> etaMin{"etaMin", -10., "Minimum eta of Particle 2 (Cascade)"};
    Configurable<float> etaMax{"etaMax", 10., "Maximum eta of Particle 2 (Cascade)"};
    Configurable<bool> useChildCuts{"useChildCuts", true, "Use cuts on the children of the Cascades additional to those of the selection of the cascade builder (for debugging purposes)"};
    Configurable<bool> useChildPIDCuts{"useChildPIDCuts", true, "Use PID cuts on the children of the Cascades additional to those of the selection of the cascade builder (for debugging purposes)"};
  } Cascade1;

  /// Partition for particle 1
  Partition<FDParticles> partitionCascade1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kCascade)) &&
                                             ((aod::femtodreamparticle::cut & Cascade1.cutBit) == Cascade1.cutBit) &&
                                             (aod::femtodreamparticle::pt > Cascade1.ptMin) &&
                                             (aod::femtodreamparticle::pt < Cascade1.ptMax) &&
                                             (aod::femtodreamparticle::eta > Cascade1.etaMin) &&
                                             (aod::femtodreamparticle::eta < Cascade1.etaMax) &&
                                             (aod::femtodreamparticle::mLambda > Cascade1.invMassMin) &&
                                             (aod::femtodreamparticle::mLambda < Cascade1.invMassMax) &&
                                             (aod::femtodreamparticle::mAntiLambda > Cascade1.invMassV0DaughMin) &&
                                             (aod::femtodreamparticle::mAntiLambda < Cascade1.invMassV0DaughMax);
  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascade, 1> cascHistoPartOne;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeV0Child, 3> posChildHistosPartOne;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeV0Child, 4> negChildHistosPartOne;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeBachelor, 8> bachChildHistosPartOne;

  /// Particle 2 (Cascade)
  struct : ConfigurableGroup {
    std::string prefix = std::string("Cascade2");
    Configurable<int> pdgCode{"pdgCode", 3334, "PDG code of particle 2 (Cascade)"};
    Configurable<femtodreamparticle::cutContainerType> cutBit{"cutBit", 32221874, "Selection bit for particle 2 (Cascade)"};
    Configurable<femtodreamparticle::cutContainerType> childPosCutBit{"childPosCutBit", 278, "Selection bit for positive child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> childPosTPCBit{"childPosTPCBit", 1024, "PID TPC bit for positive child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> childNegCutBit{"childNegCutBit", 277, "Selection bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> childNegTPCBit{"childNegTPCBit", 4096, "PID TPC bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> childBachCutBit{"childBachCutBit", 277, "Selection bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> childBachTPCBit{"childBachTPCBit", 64, "PID TPC bit for bachelor child of Cascade"};
    Configurable<float> invMassMin{"invMassMin", 1.2, "Minimum invariant mass of Particle 2 (Cascade)"};
    Configurable<float> invMassMax{"invMassMax", 1.4, "Maximum invariant mass of Particle 2 (Cascade)"};
    Configurable<float> invMassV0DaughMin{"invMassV0DaughMin", 0., "Minimum invariant mass of the V0 Daughter"};   // (???????)
    Configurable<float> invMassV0DaughMax{"invMassV0DaughMax", 999., "Maximum invariant mass of the V0 Daughter"}; // (???????)
    Configurable<float> ptMin{"ptMin", 0., "Minimum pT of Particle 2 (Cascade)"};
    Configurable<float> ptMax{"ptMax", 999., "Maximum pT of Particle 2 (Cascade)"};
    Configurable<float> etaMin{"etaMin", -10., "Minimum eta of Particle 2 (Cascade)"};
    Configurable<float> etaMax{"etaMax", 10., "Maximum eta of Particle 2 (Cascade)"};
    Configurable<bool> useChildCuts{"useChildCuts", true, "Use cuts on the children of the Cascades additional to those of the selection of the cascade builder (for debugging purposes)"};
    Configurable<bool> useChildPIDCuts{"useChildPIDCuts", true, "Use PID cuts on the children of the Cascades additional to those of the selection of the cascade builder (for debugging purposes)"};
  } Cascade2;

  /// Partition for particle 2
  Partition<FDParticles> partitionCascade2 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kCascade)) &&
                                             ((aod::femtodreamparticle::cut & Cascade2.cutBit) == Cascade2.cutBit) &&
                                             (aod::femtodreamparticle::pt > Cascade2.ptMin) &&
                                             (aod::femtodreamparticle::pt < Cascade2.ptMax) &&
                                             (aod::femtodreamparticle::eta > Cascade2.etaMin) &&
                                             (aod::femtodreamparticle::eta < Cascade2.etaMax) &&
                                             (aod::femtodreamparticle::mLambda > Cascade2.invMassMin) &&
                                             (aod::femtodreamparticle::mLambda < Cascade2.invMassMax) &&
                                             (aod::femtodreamparticle::mAntiLambda > Cascade2.invMassV0DaughMin) &&
                                             (aod::femtodreamparticle::mAntiLambda < Cascade2.invMassV0DaughMax);
  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascade, 2> cascHistoPartTwo;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeV0Child, 9> posChildHistosPartTwo;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeV0Child, 10> negChildHistosPartTwo;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeBachelor, 11> bachChildHistosPartTwo;

  /// Binning configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("Binning");
    ConfigurableAxis tempFitVarCascade{"tempFitVarCascade", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Cascade)"};
    ConfigurableAxis tempFitVarCascadeChild{"tempFitVarCascadeChild", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Cascade child)"};
    ConfigurableAxis invMass{"invMass", {200, 1.22, 1.42}, "invMass binning"};
    ConfigurableAxis pTCascade{"pTCascade", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (Cascade)"};
    ConfigurableAxis pTCascadeChild{"pTCascadeChild", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (Cascade)"};
    ConfigurableAxis pT{"pT", {20, 0.5, 4.05}, "pT binning"};
    ConfigurableAxis kstar{"kstar", {1500, 0., 6.}, "binning kstar"};
    ConfigurableAxis kT{"kT", {150, 0., 9.}, "binning kT"};
    ConfigurableAxis mT{"mT", {225, 0., 7.5}, "binning mT"};
    ConfigurableAxis multTempFit{"multTempFit", {1, 0, 1}, "multiplicity for the TempFitVar plot"};
  } Binning;
  struct : ConfigurableGroup {
    std::string prefix = std::string("Binning4D");
    ConfigurableAxis kstar{"kstar", {1500, 0., 6.}, "binning kstar for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
    ConfigurableAxis mT{"mT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
    ConfigurableAxis mult{"mult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
    ConfigurableAxis multPercentile{"multPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
  } Binning4D;

  // Mixing configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("Mixing");
    ConfigurableAxis binMult{"binMult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "bins - multiplicity"};
    ConfigurableAxis binMultPercentile{"binMultPercentile", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "bins - multiplicity percentile"};
    ConfigurableAxis binVztx{"binVztx", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "bins - z-vertex"};
    Configurable<int> depth{"depth", 5, "Number of events for mixing"};
    Configurable<int> binPolicy{"binPolicy", 0, "Binning binPolicy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both"};
  } Mixing;
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{Mixing.binVztx, Mixing.binMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{Mixing.binVztx, Mixing.binMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{Mixing.binVztx, Mixing.binMult, Mixing.binMultPercentile}, true};
  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kCascade, aod::femtodreamparticle::ParticleType::kCascade> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kCascadeV0Child, aod::femtodreamparticle::ParticleType::kCascadeV0Child> pairCloseRejectionSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kCascadeV0Child, aod::femtodreamparticle::ParticleType::kCascadeV0Child> pairCloseRejectionME;

  static constexpr uint32_t kSignPlusMask = 1 << 1;

  /// Histogram output
  HistogramRegistry registry{"Output", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext&)
  {
    // setup binnnig binPolicy for mixing
    colBinningMult = {{Mixing.binVztx, Mixing.binMult}, true};
    colBinningMultPercentile = {{Mixing.binVztx, Mixing.binMultPercentile}, true};
    colBinningMultMultPercentile = {{Mixing.binVztx, Mixing.binMult, Mixing.binMultPercentile}, true};
    eventHisto.init(&registry, Option.isMC);

    cascHistoPartOne.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTCascade, Option.dummy, Option.dummy, Binning.tempFitVarCascade, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.isMC, Cascade1.pdgCode);
    posChildHistosPartOne.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTCascadeChild, Option.dummy, Option.dummy, Binning.tempFitVarCascadeChild, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, false, 0);
    negChildHistosPartOne.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTCascadeChild, Option.dummy, Option.dummy, Binning.tempFitVarCascadeChild, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, false, 0);
    bachChildHistosPartOne.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTCascadeChild, Option.dummy, Option.dummy, Binning.tempFitVarCascadeChild, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, false, 0);

    if (!Option.sameSpecies) {
      cascHistoPartTwo.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTCascade, Option.dummy, Option.dummy, Binning.tempFitVarCascade, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Binning.invMass, Option.dummy, Option.isMC, Cascade2.pdgCode);
      posChildHistosPartTwo.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTCascadeChild, Option.dummy, Option.dummy, Binning.tempFitVarCascadeChild, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, false, 0);
      negChildHistosPartTwo.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTCascadeChild, Option.dummy, Option.dummy, Binning.tempFitVarCascadeChild, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, false, 0);
      bachChildHistosPartTwo.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTCascadeChild, Option.dummy, Option.dummy, Binning.tempFitVarCascadeChild, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, false, 0);
    }

    sameEventCont.init(&registry,
                       Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.binMult, Mixing.binMultPercentile,
                       Binning4D.kstar, Binning4D.mT, Binning4D.mult, Binning4D.multPercentile,
                       Option.isMC, Option.use4D, Option.extendedPlots,
                       Option.highkstarCut,
                       Option.smearingByOrigin);

    sameEventCont.setPDGCodes(Cascade1.pdgCode, Cascade2.pdgCode);

    mixedEventCont.init(&registry,
                        Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.binMult, Mixing.binMultPercentile,
                        Binning4D.kstar, Binning4D.mT, Binning4D.mult, Binning4D.multPercentile,
                        Option.isMC, Option.use4D, Option.extendedPlots,
                        Option.highkstarCut,
                        Option.smearingByOrigin);

    mixedEventCont.setPDGCodes(Cascade1.pdgCode, Cascade2.pdgCode);

    pairCleaner.init(&registry);
    if (Option.cprOn.value) {
      pairCloseRejectionSE.init(&registry, &registry, Option.cprDeltaPhiMax.value, Option.cprDeltaEtaMax.value, Option.cprPlotPerRadii.value, 1, Option.cprOld.value);
      pairCloseRejectionME.init(&registry, &registry, Option.cprDeltaPhiMax.value, Option.cprDeltaEtaMax.value, Option.cprPlotPerRadii.value, 2, Option.cprOld.value, 99, true);
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  template <bool isMC, typename PartitionType, typename TableTracks, typename Collision>
  void doSameEvent(PartitionType& sliceCascade1, PartitionType& sliceCascade2, TableTracks const& parts, Collision const& col)
  {
    /// Histogramming same event
    for (auto const& casc : sliceCascade1) {
      const auto& posChild1 = parts.iteratorAt(casc.index() - 3);
      const auto& negChild1 = parts.iteratorAt(casc.index() - 2);
      const auto& bachChild1 = parts.iteratorAt(casc.index() - 1);

      // check cuts on V0 children
      if (Cascade1.useChildCuts) {
        if (!(((posChild1.cut() & Cascade1.childPosCutBit) == Cascade1.childPosCutBit) &&
              ((negChild1.cut() & Cascade1.childNegCutBit) == Cascade1.childNegCutBit) &&
              ((bachChild1.cut() & Cascade1.childBachCutBit) == Cascade1.childBachCutBit))) {
          continue;
        }
      }
      if (Cascade1.useChildPIDCuts) {
        if (!(((posChild1.pidcut() & Cascade1.childPosTPCBit) == Cascade1.childPosTPCBit) &&
              ((negChild1.pidcut() & Cascade1.childNegTPCBit) == Cascade1.childNegTPCBit) &&
              ((bachChild1.pidcut() & Cascade1.childBachTPCBit) == Cascade1.childBachTPCBit))) {
          continue;
        }
      }
      cascHistoPartOne.fillQA<isMC, false>(casc, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      posChildHistosPartOne.fillQA<false, false>(posChild1, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      negChildHistosPartOne.fillQA<false, false>(negChild1, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      bachChildHistosPartOne.fillQA<false, false>(bachChild1, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
    }
    if (!Option.sameSpecies) {
      for (auto const& casc : sliceCascade2) {
        const auto& posChild2 = parts.iteratorAt(casc.index() - 3);
        const auto& negChild2 = parts.iteratorAt(casc.index() - 2);
        const auto& bachChild2 = parts.iteratorAt(casc.index() - 1);

        // check cuts on V0 children
        if (Cascade2.useChildCuts) {
          if (!(((posChild2.cut() & Cascade2.childPosCutBit) == Cascade2.childPosCutBit) &&
                ((negChild2.cut() & Cascade2.childNegCutBit) == Cascade2.childNegCutBit) &&
                ((bachChild2.cut() & Cascade2.childBachCutBit) == Cascade2.childBachCutBit))) {
            continue;
          }
        }
        if (Cascade2.useChildPIDCuts) {
          if (!(((posChild2.pidcut() & Cascade2.childPosTPCBit) == Cascade2.childPosTPCBit) &&
                ((negChild2.pidcut() & Cascade2.childNegTPCBit) == Cascade2.childNegTPCBit) &&
                ((bachChild2.pidcut() & Cascade2.childBachTPCBit) == Cascade2.childBachTPCBit))) {
            continue;
          }
        }
        cascHistoPartTwo.fillQA<isMC, false>(casc, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        posChildHistosPartTwo.fillQA<false, false>(posChild2, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        negChildHistosPartTwo.fillQA<false, false>(negChild2, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        bachChildHistosPartTwo.fillQA<false, false>(bachChild2, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }

    /// Now build particle combinations
    for (auto const& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(sliceCascade1, sliceCascade2))) {
      const auto& posChild1 = parts.iteratorAt(p1.index() - 3);
      const auto& negChild1 = parts.iteratorAt(p1.index() - 2);
      const auto& bachChild1 = parts.iteratorAt(p1.index() - 1);
      const auto& posChild2 = parts.iteratorAt(p2.index() - 3);
      const auto& negChild2 = parts.iteratorAt(p2.index() - 2);
      const auto& bachChild2 = parts.iteratorAt(p2.index() - 1);

      // cuts on Cascade children still need to be applied
      if (Cascade1.useChildCuts) {
        if (!(((posChild1.cut() & Cascade1.childPosCutBit) == Cascade1.childPosCutBit) &&
              ((negChild1.cut() & Cascade1.childNegCutBit) == Cascade1.childNegCutBit) &&
              ((bachChild1.cut() & Cascade1.childBachCutBit) == Cascade1.childBachCutBit))) {
          continue;
        }
      }
      if (Cascade1.useChildPIDCuts) {
        if (!(((posChild1.pidcut() & Cascade1.childPosTPCBit) == Cascade1.childPosTPCBit) &&
              ((negChild1.pidcut() & Cascade1.childNegTPCBit) == Cascade1.childNegTPCBit) &&
              ((bachChild1.pidcut() & Cascade1.childBachTPCBit) == Cascade1.childBachTPCBit))) {
          continue;
        }
      }

      if (Cascade2.useChildCuts) {
        if (!(((posChild2.cut() & Cascade2.childPosCutBit) == Cascade2.childPosCutBit) &&
              ((negChild2.cut() & Cascade2.childNegCutBit) == Cascade2.childNegCutBit) &&
              ((bachChild2.cut() & Cascade2.childBachCutBit) == Cascade2.childBachCutBit))) {
          continue;
        }
      }
      if (Cascade2.useChildPIDCuts) {
        if (!(((posChild2.pidcut() & Cascade2.childPosTPCBit) == Cascade2.childPosTPCBit) &&
              ((negChild2.pidcut() & Cascade2.childNegTPCBit) == Cascade2.childNegTPCBit) &&
              ((bachChild2.pidcut() & Cascade2.childBachTPCBit) == Cascade2.childBachTPCBit))) {
          continue;
        }
      }

      // SE pair set:
      sameEventCont.setPair<isMC>(p1, p2, col.multNtr(), col.multV0M(), Option.use4D, Option.extendedPlots, Option.smearingByOrigin);
    }
  }

  // process Same Event
  void processSameEvent(FilteredCollision const& col, FDParticles const& parts)
  {
    eventHisto.fillQA<false>(col);
    auto sliceCascade1 = partitionCascade1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto sliceCascade2 = partitionCascade2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<false>(sliceCascade1, sliceCascade2, parts, col);
  }
  PROCESS_SWITCH(FemtoDreamPairCascadeCascade, processSameEvent, "Enable processing same event", true);

  // Mixed events
  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent(CollisionType const& cols, PartType const& parts, PartitionType& part1, PartitionType& part2, BinningType binPolicy)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(binPolicy, Mixing.depth.value, -1, cols, cols)) {
      // make sure that tracks in same events are not mixed
      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }
      auto sliceCasc1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto sliceCasc2 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceCasc1, sliceCasc2))) {
        const auto& posChild1 = parts.iteratorAt(p1.index() - 3);
        const auto& negChild1 = parts.iteratorAt(p1.index() - 2);
        const auto& bachChild1 = parts.iteratorAt(p1.index() - 1);
        const auto& posChild2 = parts.iteratorAt(p2.index() - 3);
        const auto& negChild2 = parts.iteratorAt(p2.index() - 2);
        const auto& bachChild2 = parts.iteratorAt(p2.index() - 1);
        // check cuts on Cascade children
        if (Cascade1.useChildCuts) {
          if (!(((posChild1.cut() & Cascade1.childPosCutBit) == Cascade1.childPosCutBit) &&
                ((negChild1.cut() & Cascade1.childNegCutBit) == Cascade1.childNegCutBit) &&
                ((bachChild1.cut() & Cascade1.childBachCutBit) == Cascade1.childBachCutBit))) {
            continue;
          }
        }
        if (Cascade1.useChildPIDCuts) {
          if (!(((posChild1.pidcut() & Cascade1.childPosTPCBit) == Cascade1.childPosTPCBit) &&
                ((negChild1.pidcut() & Cascade1.childNegTPCBit) == Cascade1.childNegTPCBit) &&
                ((bachChild1.pidcut() & Cascade1.childBachTPCBit) == Cascade1.childBachTPCBit))) {
            continue;
          }
        }

        if (Cascade2.useChildCuts) {
          if (!(((posChild2.cut() & Cascade2.childPosCutBit) == Cascade2.childPosCutBit) &&
                ((negChild2.cut() & Cascade2.childNegCutBit) == Cascade2.childNegCutBit) &&
                ((bachChild2.cut() & Cascade2.childBachCutBit) == Cascade2.childBachCutBit))) {
            continue;
          }
        }
        if (Cascade2.useChildPIDCuts) {
          if (!(((posChild2.pidcut() & Cascade2.childPosTPCBit) == Cascade2.childPosTPCBit) &&
                ((negChild2.pidcut() & Cascade2.childNegTPCBit) == Cascade2.childNegTPCBit) &&
                ((bachChild2.pidcut() & Cascade2.childBachTPCBit) == Cascade2.childBachTPCBit))) {
            continue;
          }
        }

        mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), Option.use4D, Option.extendedPlots, Option.smearingByOrigin);
      }
    }
  }

  // process Mixed Event
  void processMixedEvent(FilteredCollisions const& cols, FDParticles const& parts)
  {
    switch (Mixing.binPolicy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent<false>(cols, parts, partitionCascade1, partitionCascade2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<false>(cols, parts, partitionCascade1, partitionCascade2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<false>(cols, parts, partitionCascade1, partitionCascade2, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(FemtoDreamPairCascadeCascade, processMixedEvent, "Enable processing mixed events", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoDreamPairCascadeCascade>(cfgc),
  };
  return workflow;
}
