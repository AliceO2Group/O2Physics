// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamPairTaskV0V0.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Bianca Popa, TU München, bianca.popa@tum.de

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamContainer.h"
#include "PWGCF/FemtoDream/Core/femtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include "TRandom3.h"

#include <bitset>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct femtoDreamPairTaskV0V0 {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  /// General options
  struct : ConfigurableGroup {
    std::string prefix = std::string("Option");
    Configurable<bool> IsMC{"IsMC", false, "Enable additional Histogramms in the case of runninger over Monte Carlo"};
    Configurable<bool> Use4D{"Use4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
    Configurable<bool> ExtendedPlots{"ExtendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
    Configurable<float> HighkstarCut{"HighkstarCut", -1., "Set a cut for high k*, above which the pairs are rejected. Set it to -1 to deactivate it"};
    Configurable<bool> SameSpecies{"SameSpecies", true, "Set to true if particle 1 and particle 2 are the same species"};
    Configurable<bool> MixEventWithPairs{"MixEventWithPairs", false, "Only use events that contain particle 1 and partile 2 for the event mixing"};
    Configurable<bool> RandomizePair{"RandomizePair", false, "Randomly mix particle 1 and particle 2 in case both are identical"};
    Configurable<bool> CPROn{"CPROn", true, "Close Pair Rejection"};
    Configurable<bool> CPROld{"CPROld", false, "Set to FALSE to use fixed version of CPR (for testing now, will be default soon)"};
    Configurable<bool> CPRSepMeSe{"CPRSepMESE", true, "Use seperated plots for same and mixed event for CPR plots"};
    Configurable<bool> CPRPlotPerRadii{"CPRPlotPerRadii", false, "Plot CPR per radii"};
    Configurable<float> CPRdeltaPhiMax{"CPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
    Configurable<float> CPRdeltaEtaMax{"CPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
    Configurable<bool> DCACutPtDep{"DCACutPtDep", false, "Use pt dependent dca cut"};
    Configurable<bool> SmearingByOrigin{"SmearingByOrigin", false, "Obtain the smearing matrix differential in the MC origin of particle 1 and particle 2. High memory consumption"};
    ConfigurableAxis Dummy{"Dummy", {1, 0, 1}, "Dummy axis"};
  } Option;

  /// Event selection
  struct : ConfigurableGroup {
    std::string prefix = std::string("EventSel");
    Configurable<int> MultMin{"MultMin", 0, "Minimum Multiplicity (MultNtr)"};
    Configurable<int> MultMax{"MultMax", 99999, "Maximum Multiplicity (MultNtr)"};
    Configurable<float> MultPercentileMin{"MultPercentileMin", 0, "Maximum Multiplicity Percentile"};
    Configurable<float> MultPercentileMax{"MultPercentileMax", 100, "Minimum Multiplicity Percentile"};
  } EventSel;

  Filter EventMultiplicity = aod::femtodreamcollision::multNtr >= EventSel.MultMin && aod::femtodreamcollision::multNtr <= EventSel.MultMax;
  Filter EventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= EventSel.MultPercentileMin && aod::femtodreamcollision::multV0M <= EventSel.MultPercentileMax;

  using FilteredCollisions = soa::Filtered<FDCollisions>;
  using FilteredCollision = FilteredCollisions::iterator;
  using FilteredMCCollisions = soa::Filtered<soa::Join<aod::FDCollisions, aod::FDMCCollLabels>>;
  using FilteredMCCollision = FilteredMCCollisions::iterator;

  using FilteredMaskedCollisions = soa::Filtered<soa::Join<FDCollisions, FDColMasks, FDDownSample>>;
  using FilteredMaskedCollision = FilteredMaskedCollisions::iterator;
  using FilteredMaskedMCCollisions = soa::Filtered<soa::Join<FDCollisions, aod::FDMCCollLabels, FDColMasks, FDDownSample>>;
  using FilteredMaskedMCCollision = FilteredMaskedMCCollisions::iterator;

  femtodreamcollision::BitMaskType BitMask = 0;

  /// Particle 1 (V0)
  struct : ConfigurableGroup {
    std::string prefix = std::string("V01");
    Configurable<int> PDGCode{"PDGCode", 3122, "PDG code of particle 1 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> CutBit{"CutBit", 7518, "Selection bit for particle 1 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_CutBit{"ChildPos_CutBit", 210, "Selection bit for positive child of V0"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_TPCBit{"ChildPos_TPCBit", 64, "PID TPC bit for positive child of V0"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_CutBit{"ChildNeg_CutBit", 209, "Selection bit for negative child of V0"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_TPCBit{"ChildNeg_TPCBit", 256, "PID TPC bit for negative child of V0"};

    Configurable<float> InvMassMin{"InvMassMin", 1.08, "Minimum invariant mass of Partricle 1 (particle) (V0)"};
    Configurable<float> InvMassMax{"InvMassMax", 1.15, "Maximum invariant mass of Partricle 1 (particle) (V0)"};
    Configurable<float> InvMassAntiMin{"InvMassAntiMin", 0., "Minimum invariant mass of Partricle 1 (antiparticle) (V0)"};
    Configurable<float> InvMassAntiMax{"InvMassAntiMax", 999., "Maximum invariant mass of Partricle 1 (antiparticle) (V0)"};

    Configurable<float> PtMin{"PtMin", 0., "Minimum pT of Partricle 1 (V0)"};
    Configurable<float> PtMax{"PtMax", 999., "Maximum pT of Partricle 1 (V0)"};
    Configurable<float> EtaMin{"EtaMin", -10., "Minimum eta of Partricle 1 (V0)"};
    Configurable<float> EtaMax{"EtaMax", 10., "Maximum eta of Partricle 1 (V0)"};
  } V01;

  /// Partition for particle 1
  Partition<FDParticles> PartitionV01 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                        ((aod::femtodreamparticle::cut & V01.CutBit) == V01.CutBit) &&
                                        (aod::femtodreamparticle::pt > V01.PtMin) &&
                                        (aod::femtodreamparticle::pt < V01.PtMax) &&
                                        (aod::femtodreamparticle::eta > V01.EtaMin) &&
                                        (aod::femtodreamparticle::eta < V01.EtaMax) &&
                                        (aod::femtodreamparticle::mLambda > V01.InvMassMin) &&
                                        (aod::femtodreamparticle::mLambda < V01.InvMassMax) &&
                                        (aod::femtodreamparticle::mAntiLambda > V01.InvMassAntiMin) &&
                                        (aod::femtodreamparticle::mAntiLambda < V01.InvMassAntiMax);

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 1> trackHistoPartOne;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 4> negChildHistos;

  /// Particle 2 (V0)
  struct : ConfigurableGroup {
    std::string prefix = std::string("V02");
    Configurable<int> PDGCode{"PDGCode", 3122, "PDG code of particle 2 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> CutBit{"CutBit", 7518, "Selection bit for particle 2 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_CutBit{"ChildPos_CutBit", 210, "Selection bit for positive child of V0"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_TPCBit{"ChildPos_TPCBit", 64, "PID TPC bit for positive child of V0"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_CutBit{"ChildNeg_CutBit", 209, "Selection bit for negative child of V0"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_TPCBit{"ChildNeg_TPCBit", 256, "PID TPC bit for negative child of V0"};

    Configurable<float> InvMassMin{"InvMassMin", 1.08, "Minimum invariant mass of Partricle 2 (particle) (V0)"};
    Configurable<float> InvMassMax{"InvMassMax", 1.15, "Maximum invariant mass of Partricle 2 (particle) (V0)"};
    Configurable<float> InvMassAntiMin{"InvMassAntiMin", 0., "Minimum invariant mass of Partricle 2 (antiparticle) (V0)"};
    Configurable<float> InvMassAntiMax{"InvMassAntiMax", 999., "Maximum invariant mass of Partricle 2 (antiparticle) (V0)"};

    Configurable<float> PtMin{"PtMin", 0., "Minimum pT of Partricle 2 (V0)"};
    Configurable<float> PtMax{"PtMax", 999., "Maximum pT of Partricle 2 (V0)"};
    Configurable<float> EtaMin{"EtaMin", -10., "Minimum eta of Partricle 2 (V0)"};
    Configurable<float> EtaMax{"EtaMax", 10., "Maximum eta of Partricle 2 (V0)"};
  } V02;

  /// Partition for particle 2
  Partition<FDParticles> PartitionV02 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                        ((aod::femtodreamparticle::cut & V02.CutBit) == V02.CutBit) &&
                                        (aod::femtodreamparticle::pt > V02.PtMin) &&
                                        (aod::femtodreamparticle::pt < V02.PtMax) &&
                                        (aod::femtodreamparticle::eta > V02.EtaMin) &&
                                        (aod::femtodreamparticle::eta < V02.EtaMax) &&
                                        (aod::femtodreamparticle::mLambda > V02.InvMassMin) &&
                                        (aod::femtodreamparticle::mLambda < V02.InvMassMax) &&
                                        (aod::femtodreamparticle::mAntiLambda > V02.InvMassAntiMin) &&
                                        (aod::femtodreamparticle::mAntiLambda < V02.InvMassAntiMax);

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// Binning configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("Binning");
    ConfigurableAxis TempFitVarTrack{"TempFitVarTrack", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Track)"};
    ConfigurableAxis TempFitVarV0{"TempFitVarV0", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0)"};
    ConfigurableAxis TempFitVarV0Child{"TempFitVarV0Child", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0 child)"};
    ConfigurableAxis InvMass{"InvMass", {200, 1, 1.2}, "InvMass binning"};
    ConfigurableAxis pTTrack{"pTTrack", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (Track)"};
    ConfigurableAxis pTV0{"pTV0", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0)"};
    ConfigurableAxis pTV0Child{"pTV0Child", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0)"};
    ConfigurableAxis pT{"pT", {20, 0.5, 4.05}, "pT binning"};
    ConfigurableAxis kstar{"kstar", {1500, 0., 6.}, "binning kstar"};
    ConfigurableAxis kT{"kT", {150, 0., 9.}, "binning kT"};
    ConfigurableAxis mT{"mT", {225, 0., 7.5}, "binning mT"};
    ConfigurableAxis multTempFit{"multTempFit", {1, 0, 1}, "multiplicity for the TempFitVar plot"};
  } Binning;

  struct : ConfigurableGroup {
    std::string prefix = "Binning4D";
    ConfigurableAxis kstar{"kstar", {1500, 0., 6.}, "binning kstar for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true)"};
    ConfigurableAxis mT{"mT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true)"};
    ConfigurableAxis mult{"mult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true)"};
    ConfigurableAxis multPercentile{"multPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true)"};
  } Binning4D;

  // Mixing configurables
  struct : ConfigurableGroup {
    std::string prefix = "Mixing";
    ConfigurableAxis MultMixBins{"MultMixBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
    ConfigurableAxis MultPercentileMixBins{"MultPercentileMixBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Mixing bins - multiplicity percentile"};
    ConfigurableAxis VztxMixBins{"VztxMixBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    Configurable<int> Depth{"Depth", 5, "Number of events for mixing"};
    Configurable<int> Policy{"Policy", 0, "Binning policy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both"};
  } Mixing;

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{Mixing.VztxMixBins, Mixing.MultMixBins}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{Mixing.VztxMixBins, Mixing.MultPercentileMixBins}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{Mixing.VztxMixBins, Mixing.MultMixBins, Mixing.MultPercentileMixBins}, true};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  // FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0> pairCloseRejectionSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0> pairCloseRejectionME;
  /// Histogram output
  HistogramRegistry Registry{"Output", {}, OutputObjHandlingPolicy::AnalysisObject};

  TRandom3* random;

  void init(InitContext&)
  {
    // setup columnpolicy for binning
    colBinningMult = {{Mixing.VztxMixBins, Mixing.MultMixBins}, true};
    colBinningMultPercentile = {{Mixing.VztxMixBins, Mixing.MultPercentileMixBins}, true};
    colBinningMultMultPercentile = {{Mixing.VztxMixBins, Mixing.MultMixBins, Mixing.MultPercentileMixBins}, true};

    if (Option.RandomizePair.value) {
      random = new TRandom3(0);
    }
    eventHisto.init(&Registry, Option.IsMC);
    trackHistoPartOne.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pT, Option.Dummy, Option.Dummy, Binning.TempFitVarV0, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Binning.InvMass, Option.Dummy, Option.IsMC, V01.PDGCode);
    posChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTV0Child, Option.Dummy, Option.Dummy, Binning.TempFitVarV0Child, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);
    negChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTV0Child, Option.Dummy, Option.Dummy, Binning.TempFitVarV0Child, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);

    sameEventCont.init(&Registry,
                       Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.MultMixBins, Mixing.MultPercentileMixBins,
                       Binning4D.kstar, Binning4D.mT, Binning4D.mult, Binning4D.multPercentile,
                       Option.IsMC, Option.Use4D, Option.ExtendedPlots,
                       Option.HighkstarCut,
                       Option.SmearingByOrigin);

    mixedEventCont.init(&Registry,
                        Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.MultMixBins, Mixing.MultPercentileMixBins,
                        Binning4D.kstar, Binning4D.mT, Binning4D.mult, Binning4D.multPercentile,
                        Option.IsMC, Option.Use4D, Option.ExtendedPlots,
                        Option.HighkstarCut,
                        Option.SmearingByOrigin);
    sameEventCont.setPDGCodes(V01.PDGCode, V02.PDGCode);
    mixedEventCont.setPDGCodes(V01.PDGCode, V02.PDGCode);
    // pairCleaner.init(&Registry);

    if (Option.CPROn.value) {
      pairCloseRejectionSE.init(&Registry, &Registry, Option.CPRdeltaPhiMax.value, Option.CPRdeltaEtaMax.value, Option.CPRPlotPerRadii.value, 1, Option.CPROld.value);
      pairCloseRejectionME.init(&Registry, &Registry, Option.CPRdeltaPhiMax.value, Option.CPRdeltaEtaMax.value, Option.CPRPlotPerRadii.value, 2, Option.CPROld.value);
    }
  };

  template <bool isMC, typename CollisionType>
  void fillCollision(CollisionType col)
  {
    eventHisto.fillQA<isMC>(col);
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupPartsOne partition for the first particle passed by the process function
  /// @param groupPartsTwo partition for the second particle passed by the process function
  /// @param parts femtoDreamParticles table (in case of Monte Carlo joined with FemtoDreamMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType, typename Collision>
  void doSameEvent(PartitionType SliceV01, PartitionType SliceV02, PartType parts, Collision col)
  {
    for (auto& v0 : SliceV01) {
      const auto& posChild = parts.iteratorAt(v0.index() - 2);
      const auto& negChild = parts.iteratorAt(v0.index() - 1);
      if (((posChild.cut() & V01.ChildPos_CutBit) == V01.ChildPos_CutBit) &&
          ((posChild.pidcut() & V01.ChildPos_TPCBit) == V01.ChildPos_TPCBit) &&
          ((negChild.cut() & V01.ChildNeg_CutBit) == V01.ChildNeg_CutBit) &&
          ((negChild.pidcut() & V01.ChildNeg_TPCBit) == V01.ChildNeg_TPCBit)) {
        trackHistoPartOne.fillQA<isMC, false>(v0, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        posChildHistos.fillQA<false, false>(posChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        negChildHistos.fillQA<false, false>(negChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }

    /// Now build the combinations
    float rand = 0.;
    if (Option.SameSpecies.value) {
      for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(SliceV01, SliceV02))) {
        const auto& posChild_1 = parts.iteratorAt(p1.index() - 2);
        const auto& negChild_1 = parts.iteratorAt(p1.index() - 1);
        const auto& posChild_2 = parts.iteratorAt(p2.index() - 2);
        const auto& negChild_2 = parts.iteratorAt(p2.index() - 1);
        if (((posChild_1.cut() & V01.ChildPos_CutBit) == V01.ChildPos_CutBit) &&
            ((posChild_1.pidcut() & V01.ChildPos_TPCBit) == V01.ChildPos_TPCBit) &&
            ((negChild_1.cut() & V01.ChildNeg_CutBit) == V01.ChildNeg_CutBit) &&
            ((negChild_1.pidcut() & V01.ChildNeg_TPCBit) == V01.ChildNeg_TPCBit) &&
            ((posChild_2.cut() & V02.ChildPos_CutBit) == V02.ChildPos_CutBit) &&
            ((posChild_2.pidcut() & V02.ChildPos_TPCBit) == V02.ChildPos_TPCBit) &&
            ((negChild_2.cut() & V02.ChildNeg_CutBit) == V02.ChildNeg_CutBit) &&
            ((negChild_2.pidcut() & V02.ChildNeg_TPCBit) == V02.ChildNeg_TPCBit)) {

          if (Option.CPROn.value) {
            if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
              continue;
            }
          }
          /*
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        */

          if (Option.RandomizePair.value) {
            rand = random->Rndm();
          }
          if (rand <= 0.5) {
            sameEventCont.setPair<isMC>(p1, p2, col.multNtr(), col.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.SmearingByOrigin);
          } else {
            sameEventCont.setPair<isMC>(p2, p1, col.multNtr(), col.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.SmearingByOrigin);
          }
        }
      }
    } else {
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceV01, SliceV02))) {
        if (Option.CPROn.value) {
          if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }
        /*
      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      */
        sameEventCont.setPair<isMC>(p1, p2, col.multNtr(), col.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.SmearingByOrigin);
      }
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEvent(FilteredCollision& col, o2::aod::FDParticles& parts)
  {
    fillCollision<false>(col);
    auto SliceV01 = PartitionV01->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceV02 = PartitionV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceV01.size() == 0 && SliceV02.size() == 0) {
      return;
    }
    doSameEvent<false>(SliceV01, SliceV02, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskV0V0, processSameEvent, "Enable processing same event", true);

  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent_NotMasked(CollisionType& cols, PartType& parts, PartitionType& part1, PartitionType& part2, BinningType policy)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(policy, Mixing.Depth.value, -1, cols, cols)) {
      auto SliceV01 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceV02 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      if (SliceV01.size() == 0 || SliceV02.size() == 0) {
        continue;
      }
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceV01, SliceV02))) {
        const auto& posChild_1 = parts.iteratorAt(p1.globalIndex() - 2);
        const auto& negChild_1 = parts.iteratorAt(p1.globalIndex() - 1);
        const auto& posChild_2 = parts.iteratorAt(p2.globalIndex() - 2);
        const auto& negChild_2 = parts.iteratorAt(p2.globalIndex() - 1);
        if (((posChild_1.cut() & V01.ChildPos_CutBit) == V01.ChildPos_CutBit) &&
            ((posChild_1.pidcut() & V01.ChildPos_TPCBit) == V01.ChildPos_TPCBit) &&
            ((negChild_1.cut() & V01.ChildNeg_CutBit) == V01.ChildNeg_CutBit) &&
            ((negChild_1.pidcut() & V01.ChildNeg_TPCBit) == V01.ChildNeg_TPCBit) &&
            ((posChild_2.cut() & V02.ChildPos_CutBit) == V02.ChildPos_CutBit) &&
            ((posChild_2.pidcut() & V02.ChildPos_TPCBit) == V02.ChildPos_TPCBit) &&
            ((negChild_2.cut() & V02.ChildNeg_CutBit) == V02.ChildNeg_CutBit) &&
            ((negChild_2.pidcut() & V02.ChildNeg_TPCBit) == V02.ChildNeg_TPCBit)) {

          if (Option.CPROn.value) {
            if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
              continue;
            }
          }
          mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.SmearingByOrigin);
        }
      }
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoDreamParticleTable
  void processMixedEvent(FilteredCollisions& cols, o2::aod::FDParticles& parts)
  {
    switch (Mixing.Policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent_NotMasked<false>(cols, parts, PartitionV01, PartitionV02, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_NotMasked<false>(cols, parts, PartitionV01, PartitionV02, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_NotMasked<false>(cols, parts, PartitionV01, PartitionV02, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskV0V0, processMixedEvent, "Enable processing mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskV0V0>(cfgc),
  };
  return workflow;
}
