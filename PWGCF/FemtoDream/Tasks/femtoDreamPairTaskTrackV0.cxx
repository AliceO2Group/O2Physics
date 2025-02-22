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

/// \file femtoDreamPairTaskTrackTrack.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include <sys/stat.h>
#include <cstdint>
#include <vector>
#include <string>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/Expressions.h"

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

struct femtoDreamPairTaskTrackV0 {
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

  Filter EventMultiplicity = aod::femtodreamcollision::multNtr >= EventSel.MultMin && aod::femtodreamcollision::multNtr <= EventSel.MultMax;
  Filter EventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= EventSel.MultPercentileMin && aod::femtodreamcollision::multV0M <= EventSel.MultPercentileMax;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  using FilteredMaskedCollisions = soa::Filtered<soa::Join<FDCollisions, FDColMasks, FDDownSample>>;
  using FilteredMaskedCollision = FilteredMaskedCollisions::iterator;
  using FilteredMaskedMCCollisions = soa::Filtered<soa::Join<FDCollisions, aod::FDMCCollLabels, FDColMasks, FDDownSample>>;
  using FilteredMaskedMCCollision = FilteredMaskedMCCollisions::iterator;

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

  Partition<FDMCParts> PartitionMCTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
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

  /// Particle 2 (V0)
  struct : ConfigurableGroup {
    std::string prefix = std::string("V02");
    Configurable<int> PDGCode{"PDGCode", 3122, "PDG code of particle 2 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> CutBit{"CutBit", 338, "Selection bit for particle 2 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_CutBit{"ChildPos_CutBit", 149, "Selection bit for positive child of V0"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_TPCBit{"ChildPos_TPCBit", 2, "PID TPC bit for positive child of V0"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_CutBit{"ChildNeg_CutBit", 149, "Selection bit for negative child of V0"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_TPCBit{"ChildNeg_TPCBit", 2, "PID TPC bit for negative child of V0"};

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

  Partition<FDMCParts> PartitionMCV02 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                        ((aod::femtodreamparticle::cut & V02.CutBit) == V02.CutBit) &&
                                        (aod::femtodreamparticle::pt > V02.PtMin) &&
                                        (aod::femtodreamparticle::pt < V02.PtMax) &&
                                        (aod::femtodreamparticle::eta > V02.EtaMax) &&
                                        (aod::femtodreamparticle::eta < V02.EtaMin) &&
                                        (aod::femtodreamparticle::mLambda > V02.InvMassMin) &&
                                        (aod::femtodreamparticle::mLambda < V02.InvMassMax) &&
                                        (aod::femtodreamparticle::mAntiLambda > V02.InvMassAntiMin) &&
                                        (aod::femtodreamparticle::mAntiLambda < V02.InvMassAntiMax);

  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 2> trackHistoPartTwo;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 4> negChildHistos;

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
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCloseRejectionSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCloseRejectionME;

  /// Histogram output
  HistogramRegistry Registry{"Output", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext& context)
  {
    // setup binnnig policy for mixing
    colBinningMult = {{Mixing.BinVztx, Mixing.BinMult}, true};
    colBinningMultPercentile = {{Mixing.BinVztx, Mixing.BinMultPercentile}, true};
    colBinningMultMultPercentile = {{Mixing.BinVztx, Mixing.BinMult, Mixing.BinMultPercentile}, true};

    eventHisto.init(&Registry, Option.IsMC);
    trackHistoPartOne.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTTrack, Option.Dummy, Option.Dummy, Binning.TempFitVarTrack, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.IsMC, Track1.PDGCode);
    trackHistoPartTwo.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTV0, Option.Dummy, Option.Dummy, Binning.TempFitVarV0, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Binning.InvMass, Option.Dummy, Option.IsMC, V02.PDGCode);
    posChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTV0Child, Option.Dummy, Option.Dummy, Binning.TempFitVarV0Child, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);
    negChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTV0Child, Option.Dummy, Option.Dummy, Binning.TempFitVarV0Child, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);

    sameEventCont.init(&Registry,
                       Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.BinMult, Mixing.BinMultPercentile,
                       Binning4D.kstar, Binning4D.mT, Binning4D.Mult, Binning4D.multPercentile,
                       Option.IsMC, Option.Use4D, Option.ExtendedPlots,
                       Option.HighkstarCut,
                       Option.smearingByOrigin);

    sameEventCont.setPDGCodes(Track1.PDGCode, V02.PDGCode);
    mixedEventCont.init(&Registry,
                        Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.BinMult, Mixing.BinMultPercentile,
                        Binning4D.kstar, Binning4D.mT, Binning4D.Mult, Binning4D.multPercentile,
                        Option.IsMC, Option.Use4D, Option.ExtendedPlots,
                        Option.HighkstarCut,
                        Option.smearingByOrigin);

    mixedEventCont.setPDGCodes(Track1.PDGCode, V02.PDGCode);
    pairCleaner.init(&Registry);
    if (Option.CPROn.value) {
      pairCloseRejectionSE.init(&Registry, &Registry, Option.CPRdeltaPhiMax.value, Option.CPRdeltaEtaMax.value, Option.CPRPlotPerRadii.value, 1, Option.CPROld.value);
      pairCloseRejectionME.init(&Registry, &Registry, Option.CPRdeltaPhiMax.value, Option.CPRdeltaEtaMax.value, Option.CPRPlotPerRadii.value, 2, Option.CPROld.value, 99, true);
    }

    // get bit for the collision mask
    std::bitset<8 * sizeof(femtodreamcollision::BitMaskType)> mask;
    int index = 0;
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.find("femto-dream-pair-task-track-v0") != std::string::npos) {
        if (containsNameValuePair(device.options, "Track1.CutBit", Track1.CutBit.value) &&
            containsNameValuePair(device.options, "Track1.TPCBit", Track1.TPCBit.value) &&
            containsNameValuePair(device.options, "Track1.TPCBit_Reject", Track1.TPCBit_Reject.value) &&
            containsNameValuePair(device.options, "Track1.TPCTOFBit", Track1.TPCTOFBit.value) &&
            containsNameValuePair(device.options, "Track1.PIDThres", Track1.PIDThres.value) &&
            containsNameValuePair(device.options, "Track1.PtMin", Track1.PtMin.value) &&
            containsNameValuePair(device.options, "Track1.PtMax", Track1.PtMax.value) &&
            containsNameValuePair(device.options, "Track1.EtaMin", Track1.EtaMin.value) &&
            containsNameValuePair(device.options, "Track1.EtaMax", Track1.EtaMax.value) &&
            containsNameValuePair(device.options, "Track1.TempFitVarMin", Track1.TempFitVarMin.value) &&
            containsNameValuePair(device.options, "Track1.TempFitVarMax", Track1.TempFitVarMax.value) &&
            containsNameValuePair(device.options, "V02.CutBit", V02.CutBit.value) &&
            containsNameValuePair(device.options, "V02.ChildPos_CutBit", V02.ChildPos_CutBit.value) &&
            containsNameValuePair(device.options, "V02.ChildPos_TPCBit", V02.ChildPos_TPCBit.value) &&
            containsNameValuePair(device.options, "V02.ChildNeg_CutBit", V02.ChildNeg_CutBit.value) &&
            containsNameValuePair(device.options, "V02.ChildNeg_TPCBit", V02.ChildNeg_TPCBit.value) &&
            containsNameValuePair(device.options, "V02.InvMassMin", V02.InvMassMin.value) &&
            containsNameValuePair(device.options, "V02.InvMassMax", V02.InvMassMax.value) &&
            containsNameValuePair(device.options, "V02.InvMassAntiMin", V02.InvMassAntiMin.value) &&
            containsNameValuePair(device.options, "V02.InvMassAntiMax", V02.InvMassAntiMax.value) &&
            containsNameValuePair(device.options, "V02.PtMin", V02.PtMin.value) &&
            containsNameValuePair(device.options, "V02.PtMax", V02.PtMax.value) &&
            containsNameValuePair(device.options, "V02.EtaMin", V02.EtaMin.value) &&
            containsNameValuePair(device.options, "V02.EtaMax", V02.EtaMax.value)) {
          mask.set(index);
          BitMask = static_cast<femtodreamcollision::BitMaskType>(mask.to_ulong());
          LOG(info) << "Device name matched: " << device.name;
          LOG(info) << "Bitmask for collisions: " << mask.to_string();
          break;
        } else {
          index++;
        }
      }
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  template <bool isMC, typename PartitionType, typename TableTracks, typename Collision>
  void doSameEvent(PartitionType& SliceTrk1, PartitionType& SliceV02, TableTracks const& parts, Collision const& col)
  {
    /// Histogramming same event
    for (auto const& part : SliceTrk1) {
      trackHistoPartOne.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
    }
    for (auto& v0 : SliceV02) {
      const auto& posChild = parts.iteratorAt(v0.index() - 2);
      const auto& negChild = parts.iteratorAt(v0.index() - 1);
      // This is how it is supposed to work but there seems to be an issue
      // with partitions and accessing elements in tables that have been declared
      // with an SELF_INDEX column. Under investigation. Maybe need to change
      // femtdream dataformat to take special care of v0 candidates
      // auto posChild = v0.template children_as<S>().front();
      // auto negChild = v0.template children_as<S>().back();
      // check cuts on V0 children
      if (((posChild.cut() & V02.ChildPos_CutBit) == V02.ChildPos_CutBit) &&
          ((posChild.pidcut() & V02.ChildPos_TPCBit) == V02.ChildPos_TPCBit) &&
          ((negChild.cut() & V02.ChildNeg_CutBit) == V02.ChildNeg_CutBit) &&
          ((negChild.pidcut() & V02.ChildNeg_TPCBit) == V02.ChildNeg_TPCBit)) {
        trackHistoPartTwo.fillQA<isMC, false>(v0, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        posChildHistos.fillQA<false, false>(posChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        negChildHistos.fillQA<false, false>(negChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }
    /// Now build particle combinations
    for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceV02))) {
      const auto& posChild = parts.iteratorAt(p2.index() - 2);
      const auto& negChild = parts.iteratorAt(p2.index() - 1);
      // cuts on V0 children still need to be applied
      if (((posChild.cut() & V02.ChildPos_CutBit) == V02.ChildPos_CutBit) &&
          ((posChild.pidcut() & V02.ChildPos_TPCBit) == V02.ChildPos_TPCBit) &&
          ((negChild.cut() & V02.ChildNeg_CutBit) == V02.ChildNeg_CutBit) &&
          ((negChild.pidcut() & V02.ChildNeg_TPCBit) == V02.ChildNeg_TPCBit)) {
        if (Option.CPROn.value) {
          if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        sameEventCont.setPair<isMC>(p1, p2, col.multNtr(), col.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.smearingByOrigin);
      }
    }
  }

  void processSameEventMasked(FilteredMaskedCollision const& col, FDParticles const& parts)
  {
    if ((col.bitmaskTrackOne() & BitMask) != BitMask || (col.bitmaskTrackTwo() & BitMask) != BitMask) {
      return;
    }
    eventHisto.fillQA<false>(col);
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceV02 = PartitionV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<false>(SliceTrk1, SliceV02, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEventMasked, "Enable processing same event with masks", true);

  void processSameEventMCMasked(FilteredMaskedMCCollision const& col, o2::aod::FDMCCollisions&, FDMCParts const& parts, o2::aod::FDMCParticles const&)
  {
    if ((col.bitmaskTrackOne() & BitMask) != BitMask && (col.bitmaskTrackTwo() & BitMask) != BitMask) {
      return;
    }
    eventHisto.fillQA<true>(col);
    auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceMCV02 = PartitionMCV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<true>(SliceMCTrk1, SliceMCV02, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEventMCMasked, "Enable processing same event MC with masks", false);

  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent_Masked(CollisionType const& cols, PartType const& parts, PartitionType& part1, PartitionType& part2, BinningType policy)
  {

    if (Option.MixEventWithPairs.value) {
      Partition<CollisionType> PartitionMaskedCol = ncheckbit(aod::femtodreamcollision::bitmaskTrackOne, BitMask) && ncheckbit(aod::femtodreamcollision::bitmaskTrackTwo, BitMask) && aod::femtodreamcollision::downsample == true;
      PartitionMaskedCol.bindTable(cols);
      // use *Partition.mFiltered when passing the partition to mixing object
      // there is an issue when the partition is passed directly
      // workaround for now, change back once it is fixed
      for (auto const& [collision1, collision2] : selfCombinations(policy, Mixing.Depth.value, -1, *PartitionMaskedCol.mFiltered, *PartitionMaskedCol.mFiltered)) {
        // make sure that tracks in same events are not mixed
        if (collision1.globalIndex() == collision2.globalIndex()) {
          continue;
        }
        auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto SliceV02 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
        for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceV02))) {
          const auto& posChild = parts.iteratorAt(p2.index() - 2);
          const auto& negChild = parts.iteratorAt(p2.index() - 1);
          // check cuts on V0 children
          if (((posChild.cut() & V02.ChildPos_CutBit) == V02.ChildPos_CutBit) &&
              ((posChild.pidcut() & V02.ChildPos_TPCBit) == V02.ChildPos_TPCBit) &&
              ((negChild.cut() & V02.ChildNeg_CutBit) == V02.ChildNeg_CutBit) &&
              ((negChild.pidcut() & V02.ChildNeg_TPCBit) == V02.ChildNeg_TPCBit)) {
            continue;
          }
          if (Option.CPROn.value) {
            if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
              continue;
            }
          }
          if (!pairCleaner.isCleanPair(p1, p2, parts)) {
            continue;
          }
          mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.smearingByOrigin);
        }
      }
    } else {
      Partition<CollisionType> PartitionMaskedCol1 = ncheckbit(aod::femtodreamcollision::bitmaskTrackOne, BitMask) && aod::femtodreamcollision::downsample == true;
      Partition<CollisionType> PartitionMaskedCol2 = ncheckbit(aod::femtodreamcollision::bitmaskTrackTwo, BitMask) && aod::femtodreamcollision::downsample == true;
      PartitionMaskedCol1.bindTable(cols);
      PartitionMaskedCol2.bindTable(cols);

      // use *Partition.mFiltered when passing the partition to mixing object
      // there is an issue when the partition is passed directly
      // workaround for now, change back once it is fixed
      for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockUpperIndexPolicy(policy, Mixing.Depth.value, -1, *PartitionMaskedCol1.mFiltered, *PartitionMaskedCol2.mFiltered))) {
        // make sure that tracks in same events are not mixed
        if (collision1.globalIndex() == collision2.globalIndex()) {
          continue;
        }
        auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto SliceV02 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
        for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceV02))) {
          const auto& posChild = parts.iteratorAt(p2.index() - 2);
          const auto& negChild = parts.iteratorAt(p2.index() - 1);
          // check cuts on V0 children
          if (((posChild.cut() & V02.ChildPos_CutBit) == V02.ChildPos_CutBit) &&
              ((posChild.pidcut() & V02.ChildPos_TPCBit) == V02.ChildPos_TPCBit) &&
              ((negChild.cut() & V02.ChildNeg_CutBit) == V02.ChildNeg_CutBit) &&
              ((negChild.pidcut() & V02.ChildNeg_TPCBit) == V02.ChildNeg_TPCBit)) {
            continue;
          }
          if (Option.CPROn.value) {
            if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
              continue;
            }
          }
          if (!pairCleaner.isCleanPair(p1, p2, parts)) {
            continue;
          }
          mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.smearingByOrigin);
        }
      }
    }
  }

  void processMixedEventMasked(FilteredMaskedCollisions const& cols, FDParticles const& parts)
  {
    switch (Mixing.Policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent_Masked<false>(cols, parts, PartitionTrk1, PartitionV02, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_Masked<false>(cols, parts, PartitionTrk1, PartitionV02, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_Masked<false>(cols, parts, PartitionTrk1, PartitionV02, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processMixedEventMasked, "Enable processing mixed events with masks", true);

  void processMixedEventMCMasked(FilteredMaskedMCCollisions const& cols, o2::aod::FDMCCollisions&, FDMCParts const& parts, o2::aod::FDMCParticles const&)
  {
    switch (Mixing.Policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent_Masked<true>(cols, parts, PartitionMCTrk1, PartitionMCV02, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_Masked<true>(cols, parts, PartitionMCTrk1, PartitionMCV02, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_Masked<true>(cols, parts, PartitionMCTrk1, PartitionMCV02, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processMixedEventMCMasked, "Enable processing mixed events MC with masks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskTrackV0>(cfgc),
  };
  return workflow;
}
