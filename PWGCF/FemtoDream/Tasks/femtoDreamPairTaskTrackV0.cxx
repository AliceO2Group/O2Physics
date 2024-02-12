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

#include <Framework/Expressions.h>
#include <sys/stat.h>
#include <cstdint>
#include <vector>
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

struct femtoDreamPairTaskTrackV0 {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  /// General options
  Configurable<bool> ConfOptIsMC{"ConfOptIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<bool> ConfOptUse4D{"ConfOptUse4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
  Configurable<bool> ConfOptExtendedPlots{"ConfOptExtendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
  Configurable<float> ConfOptHighkstarCut{"ConfOptHighkstarCut", -1., "Set a cut for high k*, above which the pairs are rejected. Set it to -1 to deactivate it"};
  Configurable<bool> ConfOptUseCPR{"ConfOptCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfOptCPRPlotPerRadii{"ConfOptCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfOptCPRdeltaPhiMax{"ConfOptCPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<float> ConfOptCPRdeltaEtaMax{"ConfOptCPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
  ConfigurableAxis ConfOptDummy{"ConfOptDummy", {1, 0, 1}, "Dummy axis"};

  /// Event selection
  Configurable<int> ConfEvent_minMult{"ConfEvent_minMult", 0, "Minimum Multiplicity (MultNtr)"};
  Configurable<int> ConfEvent_maxMult{"ConfEvent_maxMult", 99999, "Maximum Multiplicity (MultNtr)"};
  Configurable<float> ConfEvent_minMultPercentile{"ConfEvent_minMultPercentile", 0, "Minimum Multiplicity Percentile"};
  Configurable<float> ConfEvent_maxMultPercentile{"ConfEvent_maxMultPercentile", 100, "Maximum Multiplicity Percentile"};

  Filter EventMultiplicity = aod::femtodreamcollision::multNtr >= ConfEvent_minMult && aod::femtodreamcollision::multNtr <= ConfEvent_maxMult;
  Filter EventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= ConfEvent_minMultPercentile && aod::femtodreamcollision::multV0M <= ConfEvent_maxMultPercentile;

  using FilteredCollisions = soa::Filtered<FDCollisions>;
  using FilteredCollision = FilteredCollisions::iterator;
  using FilteredMaskedCollisions = soa::Filtered<soa::Join<FDCollisions, FDColMasks, FDDownSample>>;
  using FilteredMaskedCollision = FilteredMaskedCollisions::iterator;
  femtodreamcollision::BitMaskType BitMask = -1;

  /// Particle 1 (track)
  Configurable<int> ConfTrk1_PDGCode{"ConfTrk1_PDGCode", 2212, "PDG code of Particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk1_CutBit{"ConfTrk1_CutBit", 5542474, "Particle 1 (Track) - Selection bit from cutCulator"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk1_TPCBit{"ConfTrk1_TPCBit", 4, "PID TPC bit from cutCulator for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk1_TPCTOFBit{"ConfTrk1_TPCTOFBit", 2, "PID TPCTOF bit from cutCulator for particle 1 (Track)"};
  Configurable<float> ConfTrk1_PIDThres{"ConfTrk1_PIDThres", 0.75, "Momentum threshold for PID selection for particle 1 (Track)"};
  Configurable<float> ConfTrk1_minPt{"ConfTrk1_minPt", 0., "Minimum pT of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_maxPt{"ConfTrk1_maxPt", 999., "Maximum pT of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_minEta{"ConfTrk1_minEta", -10., "Minimum eta of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_maxEta{"ConfTrk1_maxEta", 10., "Maximum eta of partricle 1 (Track)"};

  Filter trackPtFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt > ConfTrk1_minPt, true);
  Filter trackPtFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt < ConfTrk1_maxPt, true);
  Filter trackEtaFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta > ConfTrk1_minEta, true);
  Filter trackEtaFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta < ConfTrk1_maxEta, true);

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2 (V0)
  Configurable<int> ConfV02_PDGCode{"ConfV02_PDGCode", 3122, "PDG code of particle 2 (V0)"};
  Configurable<femtodreamparticle::cutContainerType> ConfV02_CutBit{"ConfV02_CutBit", 338, "Selection bit for particle 2 (V0)"};
  Configurable<femtodreamparticle::cutContainerType> ConfV02_ChildPos_CutBit{"ConfV02_ChildPos_CutBit", 149, "Selection bit for positive child of V0"};
  Configurable<femtodreamparticle::cutContainerType> ConfV02_ChildPos_TPCBit{"ConfV02_ChildPos_TPCBit", 2, "PID TPC bit for positive child of V0"};
  Configurable<femtodreamparticle::cutContainerType> ConfV02_ChildNeg_CutBit{"ConfV02_ChildNeg_CutBit", 149, "Selection bit for negative child of V0"};
  Configurable<femtodreamparticle::cutContainerType> ConfV02_ChildNeg_TPCBit{"ConfV02_ChildNeg_TPCBit", 2, "PID TPC bit for negative child of V0"};

  Configurable<float> ConfV02_minInvMass{"ConfV02_minInvMass", 1.08, "Minimum invariant mass of Partricle 2 (particle) (V0)"};
  Configurable<float> ConfV02_maxInvMass{"ConfV02_maxInvMass", 1.15, "Maximum invariant mass of Partricle 2 (particle) (V0)"};
  Configurable<float> ConfV02_minInvMassAnti{"ConfV02_minInvMassAnti", 0., "Minimum invariant mass of Partricle 2 (antiparticle) (V0)"};
  Configurable<float> ConfV02_maxInvMassAnti{"ConfV02_maxInvMassAnti", 999., "Maximum invariant mass of Partricle 2 (antiparticle) (V0)"};

  Configurable<float> ConfV02_minPt{"ConfV02_minPt", 0., "Minimum pT of Partricle 2 (V0)"};
  Configurable<float> ConfV02_maxPt{"ConfV02_maxPt", 999., "Maximum pT of Partricle 2 (V0)"};
  Configurable<float> ConfV02_minEta{"ConfV02_minEta", -10., "Minimum eta of Partricle 2 (V0)"};
  Configurable<float> ConfV02_maxEta{"ConfV02_maxEta", 10., "Maximum eta of Partricle 2 (V0)"};

  Filter v0MassFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0), aod::femtodreamparticle::mLambda > ConfV02_minInvMass, true);
  Filter v0MassFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0), aod::femtodreamparticle::mLambda < ConfV02_maxInvMass, true);
  Filter antiv0MassFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0), aod::femtodreamparticle::mAntiLambda > ConfV02_minInvMassAnti, true);
  Filter antiv0MassFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0), aod::femtodreamparticle::mAntiLambda < ConfV02_maxInvMassAnti, true);

  Filter v0PtFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0), aod::femtodreamparticle::pt > ConfV02_minPt, true);
  Filter v0PtFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0), aod::femtodreamparticle::pt < ConfV02_maxPt, true);
  Filter v0EtaFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0), aod::femtodreamparticle::eta > ConfV02_minEta, true);
  Filter v0EtaFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0), aod::femtodreamparticle::eta < ConfV02_maxEta, true);

  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 2> trackHistoPartTwo;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 4> negChildHistos;

  using FilteredFDParticles = soa::Filtered<aod::FDParticles>;
  using FilteredFDParticle = FilteredFDParticles::iterator;

  using FilteredFDMCParts = soa::Filtered<soa::Join<aod::FDParticles, aod::FDMCLabels>>;
  using FilteredFDMCPart = FilteredFDMCParts::iterator;

  /// Partition for particle 1
  Partition<FilteredFDParticles> PartitionTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 (ncheckbit(aod::femtodreamparticle::cut, ConfTrk1_CutBit)) &&
                                                 ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk1_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCTOFBit));
  Partition<FilteredFDMCParts> PartitionMCTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 (ncheckbit(aod::femtodreamparticle::cut, ConfTrk1_CutBit)) &&
                                                 ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk1_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCTOFBit));

  /// Partition for particle 2
  Partition<FilteredFDParticles> PartitionV02 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) && ((aod::femtodreamparticle::cut & ConfV02_CutBit) == ConfV02_CutBit);
  Partition<FilteredFDMCParts> PartitionMCV02 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) && ((aod::femtodreamparticle::cut & ConfV02_CutBit) == ConfV02_CutBit);

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// Binning configurables
  ConfigurableAxis ConfBinTempFitVarTrack{"ConfBinTempFitVarTrack", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Track)"};
  ConfigurableAxis ConfBinTempFitVarV0{"ConfBinTempFitVarV0", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis ConfBinTempFitVarV0Child{"ConfBinTempFitVarV0Child", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0 child)"};
  ConfigurableAxis ConfBinInvMass{"ConfBinInvMass", {200, 1, 1.2}, "InvMass binning"};
  ConfigurableAxis ConfBinpTTrack{"ConfBinpTTrack", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (Track)"};
  ConfigurableAxis ConfBinpTV0{"ConfBinpTV0", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis ConfBinpTV0Child{"ConfBinpTV0Child", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis ConfBinpT{"ConfBinpT", {20, 0.5, 4.05}, "pT binning"};
  ConfigurableAxis ConfBinkstar{"ConfBinkstar", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfBinkT{"ConfBinkT", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfBinmT{"ConfBinmT", {225, 0., 7.5}, "binning mT"};

  ConfigurableAxis ConfBin4Dkstar{"ConfBin4Dkstar", {1500, 0., 6.}, "binning kstar for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
  ConfigurableAxis ConfBin4DmT{"ConfBin4DmT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
  ConfigurableAxis ConfBin4DMult{"ConfBin4Dmult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
  ConfigurableAxis ConfBin4DmultPercentile{"ConfBin4DmultPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};

  // Mixing configurables
  ConfigurableAxis ConfMixingBinMult{"ConfMixingBinMult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfMixingBinMultPercentile{"ConfMixingBinMultPercentile", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "Mixing bins - multiplicity percentile"};
  ConfigurableAxis ConfMixingBinVztx{"ConfMixingBinVztx", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  Configurable<int> ConfMixingDepth{"ConfMixingDepth", 5, "Number of events for mixing"};
  Configurable<int> ConfMixingPolicy{"ConfMixingBinPolicy", 0, "Binning policy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{ConfMixingBinVztx, ConfMixingBinMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{ConfMixingBinVztx, ConfMixingBinMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{ConfMixingBinVztx, ConfMixingBinMult, ConfMixingBinMultPercentile}, true};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCloseRejection;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext& context)
  {
    eventHisto.init(&qaRegistry);
    // void init(HistogramRegistry* registry,
    //    T& MomentumBins, T& tempFitVarBins, T& NsigmaTPCBins, T& NsigmaTOFBins, T& NsigmaTPCTOFBins, T& InvMassBins,
    //    bool isMC, int pdgCode, bool isDebug = false)
    trackHistoPartOne.init(&qaRegistry, ConfOptDummy, ConfOptDummy, ConfBinpTTrack, ConfOptDummy, ConfOptDummy, ConfBinTempFitVarTrack, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptIsMC, ConfTrk1_PDGCode);
    trackHistoPartTwo.init(&qaRegistry, ConfOptDummy, ConfOptDummy, ConfBinpTV0, ConfOptDummy, ConfOptDummy, ConfBinTempFitVarV0, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfBinInvMass, ConfOptIsMC, ConfV02_PDGCode);
    posChildHistos.init(&qaRegistry, ConfOptDummy, ConfOptDummy, ConfBinpTV0Child, ConfOptDummy, ConfOptDummy, ConfBinTempFitVarV0Child, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptDummy, false, 0);
    negChildHistos.init(&qaRegistry, ConfOptDummy, ConfOptDummy, ConfBinpTV0Child, ConfOptDummy, ConfOptDummy, ConfBinTempFitVarV0Child, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptDummy, ConfOptDummy, false, 0);

    sameEventCont.init(&resultRegistry,
                       ConfBinkstar, ConfBinpT, ConfBinkT, ConfBinmT, ConfMixingBinMult, ConfMixingBinMultPercentile,
                       ConfBin4Dkstar, ConfBin4DmT, ConfBin4DMult, ConfBin4DmultPercentile,
                       ConfOptIsMC, ConfOptUse4D, ConfOptExtendedPlots,
                       ConfOptHighkstarCut);

    sameEventCont.setPDGCodes(ConfTrk1_PDGCode, ConfV02_PDGCode);
    mixedEventCont.init(&resultRegistry,
                        ConfBinkstar, ConfBinpT, ConfBinkT, ConfBinmT, ConfMixingBinMult, ConfMixingBinMultPercentile,
                        ConfBin4Dkstar, ConfBin4DmT, ConfBin4DMult, ConfBin4DmultPercentile,
                        ConfOptIsMC, ConfOptUse4D, ConfOptExtendedPlots,
                        ConfOptHighkstarCut);

    mixedEventCont.setPDGCodes(ConfTrk1_PDGCode, ConfV02_PDGCode);
    pairCleaner.init(&qaRegistry);
    if (ConfOptUseCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfOptCPRdeltaPhiMax.value, ConfOptCPRdeltaEtaMax.value, ConfOptCPRPlotPerRadii.value);
    }

    // get bit for the collision mask
    std::bitset<8 * sizeof(femtodreamcollision::BitMaskType)> mask;
    int index = 0;
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.find("femto-dream-pair-task-track-v0") != std::string::npos) {
        if (containsNameValuePair(device.options, "ConfTrk1_CutBit", ConfTrk1_CutBit.value) &&
            containsNameValuePair(device.options, "ConfTrk1_TPCBit", ConfTrk1_TPCBit.value) &&
            containsNameValuePair(device.options, "ConfTrk1_TPCTOFBit", ConfTrk1_TPCTOFBit.value) &&
            containsNameValuePair(device.options, "ConfTrk1_PIDThres", ConfTrk1_PIDThres.value) &&
            containsNameValuePair(device.options, "ConfTrk1_minPt", ConfTrk1_minPt.value) &&
            containsNameValuePair(device.options, "ConfTrk1_maxPt", ConfTrk1_maxPt.value) &&
            containsNameValuePair(device.options, "ConfTrk1_minEta", ConfTrk1_minEta.value) &&
            containsNameValuePair(device.options, "ConfTrk1_maxEta", ConfTrk1_maxEta.value) &&
            containsNameValuePair(device.options, "ConfV02_CutBit", ConfV02_CutBit.value) &&
            containsNameValuePair(device.options, "ConfV02_ChildPos_CutBit", ConfV02_ChildPos_CutBit.value) &&
            containsNameValuePair(device.options, "ConfV02_ChildPos_TPCBit", ConfV02_ChildPos_TPCBit.value) &&
            containsNameValuePair(device.options, "ConfV02_ChildNeg_CutBit", ConfV02_ChildNeg_CutBit.value) &&
            containsNameValuePair(device.options, "ConfV02_ChildNeg_TPCBit", ConfV02_ChildNeg_TPCBit.value) &&
            containsNameValuePair(device.options, "ConfV02_minInvMass", ConfV02_minInvMass.value) &&
            containsNameValuePair(device.options, "ConfV02_maxInvMass", ConfV02_maxInvMass.value) &&
            containsNameValuePair(device.options, "ConfV02_minInvMassAnti", ConfV02_minInvMassAnti.value) &&
            containsNameValuePair(device.options, "ConfV02_maxInvMassAnti", ConfV02_maxInvMassAnti.value) &&
            containsNameValuePair(device.options, "ConfV02_minPt", ConfV02_minPt.value) &&
            containsNameValuePair(device.options, "ConfV02_maxPt", ConfV02_maxPt.value) &&
            containsNameValuePair(device.options, "ConfV02_minEta", ConfV02_minEta.value) &&
            containsNameValuePair(device.options, "ConfV02_maxEta", ConfV02_maxEta.value)) {
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
    if ((doprocessSameEvent && doprocessSameEventMasked) ||
        (doprocessMixedEvent && doprocessMixedEventMasked) ||
        (doprocessSameEventMC && doprocessSameEventMCMasked) ||
        (doprocessMixedEventMC && doprocessMixedEventMCMasked)) {
      LOG(fatal) << "Normal and masked processing cannot be activated simultaneously!";
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  template <bool isMC, typename PartitionType, typename TableTracks, typename Collision>
  void doSameEvent(PartitionType& SliceTrk1, PartitionType& SliceV02, TableTracks& parts, Collision col)
  {
    /// Histogramming same event
    for (auto& part : SliceTrk1) {
      trackHistoPartOne.fillQA<isMC, false>(part, col, aod::femtodreamparticle::kPt);
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
      if (((posChild.cut() & ConfV02_ChildPos_CutBit) == ConfV02_ChildPos_CutBit) &&
          ((posChild.pidcut() & ConfV02_ChildPos_TPCBit) == ConfV02_ChildPos_TPCBit) &&
          ((negChild.cut() & ConfV02_ChildNeg_CutBit) == ConfV02_ChildNeg_CutBit) &&
          ((negChild.pidcut() & ConfV02_ChildNeg_TPCBit) == ConfV02_ChildNeg_TPCBit)) {
        trackHistoPartTwo.fillQA<isMC, false>(v0, col, aod::femtodreamparticle::kPt);
        posChildHistos.fillQA<false, false>(posChild, col, aod::femtodreamparticle::kPt);
        negChildHistos.fillQA<false, false>(negChild, col, aod::femtodreamparticle::kPt);
      }
    }
    /// Now build particle combinations
    for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceV02))) {
      const auto& posChild = parts.iteratorAt(p2.index() - 2);
      const auto& negChild = parts.iteratorAt(p2.index() - 1);
      // cuts on V0 children still need to be applied
      if (((posChild.cut() & ConfV02_ChildPos_CutBit) == ConfV02_ChildPos_CutBit) &&
          ((posChild.pidcut() & ConfV02_ChildPos_TPCBit) == ConfV02_ChildPos_TPCBit) &&
          ((negChild.cut() & ConfV02_ChildNeg_CutBit) == ConfV02_ChildNeg_CutBit) &&
          ((negChild.pidcut() & ConfV02_ChildNeg_TPCBit) == ConfV02_ChildNeg_TPCBit)) {
        if (ConfOptUseCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        sameEventCont.setPair<isMC>(p1, p2, col.multNtr(), col.multV0M(), ConfOptUse4D, ConfOptExtendedPlots);
      }
    }
  }

  void processSameEvent(FilteredCollision const& col, FilteredFDParticles const& parts)
  {
    eventHisto.fillQA(col);
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceV02 = PartitionV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceTrk1.size() == 0 && SliceV02.size() == 0) {
      return;
    }
    doSameEvent<false>(SliceTrk1, SliceV02, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEvent, "Enable processing same event", true);

  void processSameEventMasked(FilteredMaskedCollision const& col, FilteredFDParticles const& parts)
  {
    if ((col.bitmaskTrackOne() & BitMask) != BitMask && (col.bitmaskTrackTwo() & BitMask) != BitMask) {
      return;
    }
    eventHisto.fillQA(col);
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceV02 = PartitionV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<false>(SliceTrk1, SliceV02, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEventMasked, "Enable processing same event with masks", false);

  void processSameEventMC(FilteredCollision& col, FilteredFDMCParts& parts, o2::aod::FDMCParticles&)
  {
    eventHisto.fillQA(col);
    auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceMCV02 = PartitionMCV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceMCTrk1.size() == 0 && SliceMCV02.size() == 0) {
      return;
    }
    doSameEvent<true>(SliceMCTrk1, SliceMCV02, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEventMC, "Enable processing same event MC", false);

  void processSameEventMCMasked(FilteredMaskedCollision& col, FilteredFDMCParts& parts, o2::aod::FDMCParticles&)
  {
    if ((col.bitmaskTrackOne() & BitMask) != BitMask && (col.bitmaskTrackTwo() & BitMask) != BitMask) {
      return;
    }
    eventHisto.fillQA(col);
    auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceMCV02 = PartitionMCV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<true>(SliceMCTrk1, SliceMCV02, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEventMCMasked, "Enable processing same event MC with masks", false);

  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent_NotMasked(CollisionType& cols, PartType& parts, PartitionType& part1, PartitionType& part2, BinningType policy)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(policy, ConfMixingDepth.value, -1, cols, cols)) {
      auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceV02 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      if (SliceTrk1.size() == 0 || SliceV02.size() == 0) {
        continue;
      }
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceV02))) {
        const auto& posChild = parts.iteratorAt(p2.index() - 2);
        const auto& negChild = parts.iteratorAt(p2.index() - 1);
        // check cuts on V0 children
        if (((posChild.cut() & ConfV02_ChildPos_CutBit) == ConfV02_ChildPos_CutBit) &&
            ((posChild.pidcut() & ConfV02_ChildPos_TPCBit) == ConfV02_ChildPos_TPCBit) &&
            ((negChild.cut() & ConfV02_ChildNeg_CutBit) == ConfV02_ChildNeg_CutBit) &&
            ((negChild.pidcut() & ConfV02_ChildNeg_TPCBit) == ConfV02_ChildNeg_TPCBit)) {
          continue;
        }
        if (ConfOptUseCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, collision1.magField())) {
            continue;
          }
        }
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), ConfOptUse4D, ConfOptExtendedPlots);
      }
    }
  }

  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent_Masked(CollisionType& cols, PartType& parts, PartitionType& part1, PartitionType& part2, BinningType policy)
  {
    Partition<CollisionType> PartitionMaskedCol1 = (aod::femtodreamcollision::bitmaskTrackOne & BitMask) == BitMask && aod::femtodreamcollision::downsample == true;
    Partition<CollisionType> PartitionMaskedCol2 = (aod::femtodreamcollision::bitmaskTrackTwo & BitMask) == BitMask && aod::femtodreamcollision::downsample == true;
    PartitionMaskedCol1.bindTable(cols);
    PartitionMaskedCol2.bindTable(cols);
    for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockUpperIndexPolicy(policy, ConfMixingDepth.value, -1, PartitionMaskedCol1, PartitionMaskedCol2))) {
      auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceV02 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceV02))) {
        const auto& posChild = parts.iteratorAt(p2.index() - 2);
        const auto& negChild = parts.iteratorAt(p2.index() - 1);
        // check cuts on V0 children
        if (((posChild.cut() & ConfV02_ChildPos_CutBit) == ConfV02_ChildPos_CutBit) &&
            ((posChild.pidcut() & ConfV02_ChildPos_TPCBit) == ConfV02_ChildPos_TPCBit) &&
            ((negChild.cut() & ConfV02_ChildNeg_CutBit) == ConfV02_ChildNeg_CutBit) &&
            ((negChild.pidcut() & ConfV02_ChildNeg_TPCBit) == ConfV02_ChildNeg_TPCBit)) {
          continue;
        }
        if (ConfOptUseCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, collision1.magField())) {
            continue;
          }
        }
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), ConfOptUse4D, ConfOptExtendedPlots);
      }
    }
  }

  void processMixedEvent(FilteredCollisions& cols, FilteredFDParticles& parts)
  {
    switch (ConfMixingPolicy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent_NotMasked<false>(cols, parts, PartitionTrk1, PartitionV02, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_NotMasked<false>(cols, parts, PartitionTrk1, PartitionV02, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_NotMasked<false>(cols, parts, PartitionTrk1, PartitionV02, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processMixedEvent, "Enable processing mixed events", true);

  void processMixedEventMasked(FilteredMaskedCollisions& cols, FilteredFDParticles& parts)
  {
    switch (ConfMixingPolicy.value) {
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
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processMixedEventMasked, "Enable processing mixed events with masks", false);

  void processMixedEventMC(FilteredCollisions& cols, FilteredFDMCParts& parts, o2::aod::FDMCParticles&)
  {
    switch (ConfMixingPolicy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent_NotMasked<true>(cols, parts, PartitionMCTrk1, PartitionMCV02, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_NotMasked<true>(cols, parts, PartitionMCTrk1, PartitionMCV02, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_NotMasked<true>(cols, parts, PartitionMCTrk1, PartitionMCV02, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processMixedEventMC, "Enable processing mixed events MC", false);

  void processMixedEventMCMasked(FilteredMaskedCollisions& cols, FilteredFDMCParts& parts, o2::aod::FDMCParticles&)
  {
    switch (ConfMixingPolicy.value) {
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
