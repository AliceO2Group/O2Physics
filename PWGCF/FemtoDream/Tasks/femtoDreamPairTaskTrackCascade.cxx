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
/// \file femtoDreamPairTaskTrackTrack.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
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
struct femtoDreamPairTaskTrackCascade {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  /// General options
  struct : ConfigurableGroup {
    std::string prefix = std::string("Option");
    Configurable<bool> IsMC{"IsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<bool> Use4D{"Use4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
    Configurable<bool> ExtendedPlots{"ExtendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
    Configurable<float> HighkstarCut{"HighkstarCut", -1., "Set a cut for high k*, above which the pairs are rejected. Set it to -1 to deactivate it"};
    Configurable<bool> CPROn{"CPROn", true, "Close Pair Rejection"};
    Configurable<bool> CPROld{"CPROld", false, "Set to FALSE to use fixed version of CPR (for testing now, will be default soon)"};
    Configurable<bool> CPRPlotPerRadii{"CPRPlotPerRadii", false, "Plot CPR per radii"};
    Configurable<float> CPRdeltaPhiMax{"CPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
    Configurable<float> CPRdeltaEtaMax{"CPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
    Configurable<bool> DCACutPtDep{"DCACutPtDep", false, "Use pt dependent dca cut"};
    Configurable<bool> MixEventWithPairs{"MixEventWithPairs", false, "Only use events that contain particle 1 and partile 2 for the event mixing"};
    Configurable<bool> smearingByOrigin{"smearingByOrigin", false, "Obtain the smearing matrix differential in the MC origin of particle 1 and particle 2. High memory consumption. Use with care!"};
    ConfigurableAxis Dummy{"Dummy", {1, 0, 1}, "Dummy axis"};
  } Option;
  /// Event selection
  struct : ConfigurableGroup {
    std::string prefix = std::string("EventSel");
    Configurable<int> MultMin{"MultMin", 0, "Minimum Multiplicity (MultNtr)"};
    Configurable<int> MultMax{"MultMax", 99999, "Maximum Multiplicity (MultNtr)"};
    Configurable<float> MultPercentileMin{"MultPercentileMin", 0, "Minimum Multiplicity Percentile"};
    Configurable<float> MultPercentileMax{"MultPercentileMax", 100, "Maximum Multiplicity Percentile"};
  } EventSel;
  // Filter EventMultiplicity = aod::femtodreamcollision::multNtr >= EventSel.MultMin && aod::femtodreamcollision::multNtr <= EventSel.MultMax;
  // Filter EventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= EventSel.MultPercentileMin && aod::femtodreamcollision::multV0M <= EventSel.MultPercentileMax;
  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;
  // using FilteredCollisions = soa::Filtered<FDCollisions>;
  using FilteredCollisions = FDCollisions;
  using FilteredCollision = FilteredCollisions::iterator;
  using FDMCParts = soa::Join<aod::FDParticles, aod::FDMCLabels>;
  using FDMCPart = FDMCParts::iterator;
  femtodreamcollision::BitMaskType BitMask = 1;
  /// Particle 1 (track)
  struct : ConfigurableGroup {
    std::string prefix = std::string("Track1");
    Configurable<int> PDGCode{"PDGCode", 2212, "PDG code of Particle 1 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> CutBit{"CutBit", 5542474, "Particle 1 (Track) - Selection bit from cutCulator"};
    Configurable<femtodreamparticle::cutContainerType> TPCBit{"TPCBit", 4, "PID TPC bit from cutCulator for particle 1 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> TPCBit_Reject{"TPCBit_Reject", 0, "Reject PID TPC bit from cutCulator for particle 1 (Track). Set to 0 to turn off"};
    Configurable<femtodreamparticle::cutContainerType> TPCTOFBit{"TPCTOFBit", 2, "PID TPCTOF bit from cutCulator for particle 1 (Track)"};
    Configurable<float> PIDThres{"PIDThres", 0.75, "Momentum threshold for PID selection for particle 1 (Track)"};
    Configurable<float> PtMin{"PtMin", 0., "Minimum pT of partricle 1 (Track)"};
    Configurable<float> PtMax{"PtMax", 999., "Maximum pT of partricle 1 (Track)"};
    Configurable<float> EtaMin{"EtaMin", -10., "Minimum eta of partricle 1 (Track)"};
    Configurable<float> EtaMax{"EtaMax", 10., "Maximum eta of partricle 1 (Track)"};
    Configurable<float> TempFitVarMin{"TempFitVarMin", -10., "Minimum DCAxy of partricle 1 (Track)"};
    Configurable<float> TempFitVarMax{"TempFitVarMax", 10., "Maximum DCAxy of partricle 1 (Track)"};
  } Track1;
  /// Partition for particle 1
  Partition<FDParticles> PartitionTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                         (ncheckbit(aod::femtodreamparticle::cut, Track1.CutBit)) &&
                                         ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= Track1.PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, Track1.TPCBit) && ((aod::femtodreamparticle::pidcut & Track1.TPCBit_Reject) == 0u), ncheckbit(aod::femtodreamparticle::pidcut, Track1.TPCTOFBit)) &&
                                         (aod::femtodreamparticle::pt > Track1.PtMin) &&
                                         (aod::femtodreamparticle::pt < Track1.PtMax) &&
                                         (aod::femtodreamparticle::eta > Track1.EtaMin) &&
                                         (aod::femtodreamparticle::eta < Track1.EtaMax) &&
                                         ifnode(Option.DCACutPtDep, (nabs(aod::femtodreamparticle::tempFitVar) <= 0.0105f + (0.035f / npow(aod::femtodreamparticle::pt, 1.1f))),
                                                ((aod::femtodreamparticle::tempFitVar >= Track1.TempFitVarMin) &&
                                                 (aod::femtodreamparticle::tempFitVar <= Track1.TempFitVarMax)));
  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;
  /// Particle 2 (Cascade)
  struct : ConfigurableGroup {
    std::string prefix = std::string("Cascade2");
    Configurable<int> PDGCode{"PDGCode", 3312, "PDG code of particle 2 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> CutBit{"CutBit", 32221874, "Selection bit for particle 2 (Cascade)"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_CutBit{"ChildPos_CutBit", 278, "Selection bit for positive child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_TPCBit{"ChildPos_TPCBit", 1024, "PID TPC bit for positive child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_CutBit{"ChildNeg_CutBit", 277, "Selection bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_TPCBit{"ChildNeg_TPCBit", 4096, "PID TPC bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildBach_CutBit{"ChildBach_CutBit", 277, "Selection bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildBach_TPCBit{"ChildBach_TPCBit", 64, "PID TPC bit for bachelor child of Cascade"};
    Configurable<float> InvMassMin{"InvMassMin", 1.2, "Minimum invariant mass of Partricle 2 (particle) (Cascade)"};
    Configurable<float> InvMassMax{"InvMassMax", 1.4, "Maximum invariant mass of Partricle 2 (particle) (Cascade)"};
    Configurable<float> InvMassV0DaughMin{"InvMassV0DaugMin", 0., "Minimum invariant mass of the V0 Daughter"};
    Configurable<float> InvMassV0DaughMax{"InvMassV0DaugMax", 999., "Maximum invariant mass of the V0 Daughter"};
    Configurable<float> PtMin{"PtMin", 0., "Minimum pT of Partricle 2 (V0)"};
    Configurable<float> PtMax{"PtMax", 999., "Maximum pT of Partricle 2 (V0)"};
    Configurable<float> EtaMin{"EtaMin", -10., "Minimum eta of Partricle 2 (V0)"};
    Configurable<float> EtaMax{"EtaMax", 10., "Maximum eta of Partricle 2 (V0)"};
    Configurable<bool> UseChildCuts{"UseChildCuts", true, "Use cuts on the children of the Cascades additional to those of the selection of the cascade builder (for debugging purposes)"};
    Configurable<bool> UseChildPIDCuts{"UseChildPIDCuts", true, "Use PID cuts on the children of the Cascades additional to those of the selection of the cascade builder (for debugging purposes)"};
  } Cascade2;
  /// Partition for particle 2
  Partition<FDParticles> PartitionCascade2 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kCascade)) &&
                                             ((aod::femtodreamparticle::cut & Cascade2.CutBit) == Cascade2.CutBit) &&
                                             (aod::femtodreamparticle::pt > Cascade2.PtMin) &&
                                             (aod::femtodreamparticle::pt < Cascade2.PtMax) &&
                                             (aod::femtodreamparticle::eta > Cascade2.EtaMin) &&
                                             (aod::femtodreamparticle::eta < Cascade2.EtaMax) &&
                                             (aod::femtodreamparticle::mLambda > Cascade2.InvMassMin) &&
                                             (aod::femtodreamparticle::mLambda < Cascade2.InvMassMax) &&
                                             (aod::femtodreamparticle::mAntiLambda > Cascade2.InvMassV0DaughMin) &&
                                             (aod::femtodreamparticle::mAntiLambda < Cascade2.InvMassV0DaughMax);
  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascade, 2> trackHistoPartTwo;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeV0Child, 3> posChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeV0Child, 4> negChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeBachelor, 8> bachChildHistos;
  /// Binning configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("Binning");
    ConfigurableAxis TempFitVarTrack{"TempFitVarTrack", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Track)"};
    ConfigurableAxis TempFitVarCascade{"TempFitVarCascade", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Cascade)"};
    ConfigurableAxis TempFitVarCascadeChild{"TempFitVarCascadeChild", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Cascade child)"};
    ConfigurableAxis InvMass{"InvMass", {200, 1.22, 1.42}, "InvMass binning"};
    ConfigurableAxis pTTrack{"pTTrack", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (Track)"};
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
    ConfigurableAxis Mult{"mult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
    ConfigurableAxis multPercentile{"multPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
  } Binning4D;
  // Mixing configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("Mixing");
    ConfigurableAxis BinMult{"BinMult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "bins - multiplicity"};
    ConfigurableAxis BinMultPercentile{"BinMultPercentile", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "bins - multiplicity percentile"};
    ConfigurableAxis BinVztx{"BinVztx", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "bins - z-vertex"};
    Configurable<int> Depth{"Depth", 5, "Number of events for mixing"};
    Configurable<int> Policy{"BinPolicy", 0, "Binning policy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both"};
  } Mixing;
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{Mixing.BinVztx, Mixing.BinMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{Mixing.BinVztx, Mixing.BinMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{Mixing.BinVztx, Mixing.BinMult, Mixing.BinMultPercentile}, true};
  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCascade> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCascadeV0Child> pairCloseRejectionSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCascadeV0Child> pairCloseRejectionME;

  static constexpr uint32_t kSignPlusMask = 1 << 1;

  /// Histogram output
  HistogramRegistry Registry{"Output", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext&)
  {
    // setup binnnig policy for mixing
    colBinningMult = {{Mixing.BinVztx, Mixing.BinMult}, true};
    colBinningMultPercentile = {{Mixing.BinVztx, Mixing.BinMultPercentile}, true};
    colBinningMultMultPercentile = {{Mixing.BinVztx, Mixing.BinMult, Mixing.BinMultPercentile}, true};
    eventHisto.init(&Registry, Option.IsMC);
    trackHistoPartOne.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTTrack, Option.Dummy, Option.Dummy, Binning.TempFitVarTrack, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.IsMC, Track1.PDGCode);
    trackHistoPartTwo.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTCascade, Option.Dummy, Option.Dummy, Binning.TempFitVarCascade, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Binning.InvMass, Option.Dummy, Option.IsMC, Cascade2.PDGCode);
    posChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTCascadeChild, Option.Dummy, Option.Dummy, Binning.TempFitVarCascadeChild, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);
    negChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTCascadeChild, Option.Dummy, Option.Dummy, Binning.TempFitVarCascadeChild, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);
    bachChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTCascadeChild, Option.Dummy, Option.Dummy, Binning.TempFitVarCascadeChild, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);
    sameEventCont.init(&Registry,
                       Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.BinMult, Mixing.BinMultPercentile,
                       Binning4D.kstar, Binning4D.mT, Binning4D.Mult, Binning4D.multPercentile,
                       Option.IsMC, Option.Use4D, Option.ExtendedPlots,
                       Option.HighkstarCut,
                       Option.smearingByOrigin);

    sameEventCont.setPDGCodes(Track1.PDGCode, Cascade2.PDGCode);

    mixedEventCont.init(&Registry,
                        Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.BinMult, Mixing.BinMultPercentile,
                        Binning4D.kstar, Binning4D.mT, Binning4D.Mult, Binning4D.multPercentile,
                        Option.IsMC, Option.Use4D, Option.ExtendedPlots,
                        Option.HighkstarCut,
                        Option.smearingByOrigin);

    mixedEventCont.setPDGCodes(Track1.PDGCode, Cascade2.PDGCode);

    pairCleaner.init(&Registry);
    if (Option.CPROn.value) {
      pairCloseRejectionSE.init(&Registry, &Registry, Option.CPRdeltaPhiMax.value, Option.CPRdeltaEtaMax.value, Option.CPRPlotPerRadii.value, 1, Option.CPROld.value);
      pairCloseRejectionME.init(&Registry, &Registry, Option.CPRdeltaPhiMax.value, Option.CPRdeltaEtaMax.value, Option.CPRPlotPerRadii.value, 2, Option.CPROld.value, 99, true);
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  template <bool isMC, typename PartitionType, typename TableTracks, typename Collision>
  void doSameEvent(PartitionType& SliceTrk1, PartitionType& SliceCascade2, TableTracks const& parts, Collision const& col)
  {
    /// Histogramming same event
    for (auto const& part : SliceTrk1) {
      trackHistoPartOne.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
    }
    for (auto& casc : SliceCascade2) {
      const auto& posChild = parts.iteratorAt(casc.index() - 3);
      const auto& negChild = parts.iteratorAt(casc.index() - 2);
      const auto& bachChild = parts.iteratorAt(casc.index() - 1);
      // This is how it is supposed to work but there seems to be an issue
      // with partitions and accessing elements in tables that have been declared
      // with an SELF_INDEX column. Under investigation. Maybe need to change
      // femtdream dataformat to take special care of v0 candidates
      // auto posChild = v0.template children_as<S>().front();
      // auto negChild = v0.template children_as<S>().back();
      // check cuts on V0 children
      if (Cascade2.UseChildCuts) {
        if (!(((posChild.cut() & Cascade2.ChildPos_CutBit) == Cascade2.ChildPos_CutBit) &&
              ((negChild.cut() & Cascade2.ChildNeg_CutBit) == Cascade2.ChildNeg_CutBit) &&
              ((bachChild.cut() & Cascade2.ChildBach_CutBit) == Cascade2.ChildBach_CutBit))) {
          continue;
        }
      }
      if (Cascade2.UseChildPIDCuts) {
        if (!(((posChild.pidcut() & Cascade2.ChildPos_TPCBit) == Cascade2.ChildPos_TPCBit) &&
              ((negChild.pidcut() & Cascade2.ChildNeg_TPCBit) == Cascade2.ChildNeg_TPCBit) &&
              ((bachChild.pidcut() & Cascade2.ChildBach_TPCBit) == Cascade2.ChildBach_TPCBit))) {
          continue;
        }
      }
      trackHistoPartTwo.fillQA<isMC, false>(casc, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      posChildHistos.fillQA<false, false>(posChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      negChildHistos.fillQA<false, false>(negChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      bachChildHistos.fillQA<false, false>(bachChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
    }
    /// Now build particle combinations
    for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceCascade2))) {
      const auto& posChild = parts.iteratorAt(p2.index() - 3);
      const auto& negChild = parts.iteratorAt(p2.index() - 2);
      const auto& bachChild = parts.iteratorAt(p2.index() - 1);

      // cuts on Cascade children still need to be applied
      if (Cascade2.UseChildCuts) {
        if (!(((posChild.cut() & Cascade2.ChildPos_CutBit) == Cascade2.ChildPos_CutBit) &&
              ((negChild.cut() & Cascade2.ChildNeg_CutBit) == Cascade2.ChildNeg_CutBit) &&
              ((bachChild.cut() & Cascade2.ChildBach_CutBit) == Cascade2.ChildBach_CutBit))) {
          continue;
        }
      }
      if (Cascade2.UseChildPIDCuts) {
        if (!(((posChild.pidcut() & Cascade2.ChildPos_TPCBit) == Cascade2.ChildPos_TPCBit) &&
              ((negChild.pidcut() & Cascade2.ChildNeg_TPCBit) == Cascade2.ChildNeg_TPCBit) &&
              ((bachChild.pidcut() & Cascade2.ChildBach_TPCBit) == Cascade2.ChildBach_TPCBit))) {
          continue;
        }
      }
      if (Option.CPROn.value) {
        if ((p1.cut() & kSignPlusMask) == kSignPlusMask) {
          if (pairCloseRejectionSE.isClosePair(p1, posChild, parts, col.magField())) {
            continue;
          }
        } else {
          if (pairCloseRejectionSE.isClosePair(p1, posChild, parts, col.magField())) {
            continue;
          }
        }
      }

      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      sameEventCont.setPair<isMC>(p1, p2, col.multNtr(), col.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.smearingByOrigin);
    }
  }
  void processSameEvent(FilteredCollision const& col, FDParticles const& parts)
  {
    // if ((col.bitmaskTrackOne() & BitMask) != BitMask || (col.bitmaskTrackTwo() & BitMask) != BitMask) {
    //   return;
    // }
    eventHisto.fillQA<false>(col);
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceCascade2 = PartitionCascade2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<false>(SliceTrk1, SliceCascade2, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackCascade, processSameEvent, "Enable processing same event", true);

  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent(CollisionType const& cols, PartType const& parts, PartitionType& part1, PartitionType& part2, BinningType policy)
  {
    // Partition<CollisionType> PartitionMaskedCol = ncheckbit(aod::femtodreamcollision::bitmaskTrackOne, BitMask) && ncheckbit(aod::femtodreamcollision::bitmaskTrackTwo, BitMask);// && aod::femtodreamcollision::downsample == true;
    // PartitionMaskedCol.bindTable(cols);

    // use *Partition.mFiltered when passing the partition to mixing object
    // there is an issue when the partition is passed directly
    // workaround for now, change back once it is fixed
    for (auto const& [collision1, collision2] : soa::selfCombinations(policy, Mixing.Depth.value, -1, cols, cols)) {
      // make sure that tracks in same events are not mixed
      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }
      auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceCasc2 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceCasc2))) {
        const auto& posChild = parts.iteratorAt(p2.globalIndex() - 3);
        const auto& negChild = parts.iteratorAt(p2.globalIndex() - 2);
        const auto& bachChild = parts.iteratorAt(p2.globalIndex() - 1);
        // check cuts on Cascade children
        if (Cascade2.UseChildCuts) {
          if (!(((posChild.cut() & Cascade2.ChildPos_CutBit) == Cascade2.ChildPos_CutBit) &&
                ((negChild.cut() & Cascade2.ChildNeg_CutBit) == Cascade2.ChildNeg_CutBit) &&
                ((bachChild.cut() & Cascade2.ChildBach_CutBit) == Cascade2.ChildBach_CutBit))) {
            continue;
          }
        }
        if (Cascade2.UseChildPIDCuts) {
          if (!(((posChild.pidcut() & Cascade2.ChildPos_TPCBit) == Cascade2.ChildPos_TPCBit) &&
                ((negChild.pidcut() & Cascade2.ChildNeg_TPCBit) == Cascade2.ChildNeg_TPCBit) &&
                ((bachChild.pidcut() & Cascade2.ChildBach_TPCBit) == Cascade2.ChildBach_TPCBit))) {
            continue;
          }
        }

        if (Option.CPROn.value) {
          if ((p1.cut() & kSignPlusMask) == kSignPlusMask) {
            if (pairCloseRejectionME.isClosePair(p1, posChild, parts, collision1.magField())) {
              continue;
            }
          } else {
            if (pairCloseRejectionME.isClosePair(p1, negChild, parts, collision1.magField())) {
              continue;
            }
          }
        }
        // Pair cleaner not needed in the mixing
        // if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        //  continue;
        //}

        mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.smearingByOrigin);
      }
    }
  }

  void processMixedEvent(FilteredCollisions const& cols, FDParticles const& parts)
  {
    switch (Mixing.Policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent<false>(cols, parts, PartitionTrk1, PartitionCascade2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<false>(cols, parts, PartitionTrk1, PartitionCascade2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<false>(cols, parts, PartitionTrk1, PartitionCascade2, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackCascade, processMixedEvent, "Enable processing mixed events", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskTrackCascade>(cfgc),
  };
  return workflow;
}
