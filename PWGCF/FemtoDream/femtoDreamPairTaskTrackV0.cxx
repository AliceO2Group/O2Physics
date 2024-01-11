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
#include "FemtoDreamParticleHisto.h"
#include "FemtoDreamEventHisto.h"
#include "FemtoDreamPairCleaner.h"
#include "FemtoDreamContainer.h"
#include "FemtoDreamDetaDphiStar.h"
#include "FemtoUtils.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct femtoDreamPairTaskTrackV0 {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  using MaskedCollisions = soa::Join<FDCollisions, FDColMasks>;
  using MaskedCollision = MaskedCollisions::iterator;
  femtodreamcollision::BitMaskType MaskBit = -1;

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

  ConfigurableAxis ConfTrkTempFitVarBins{"ConfTrkTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTrkTempFitVarpTBins{"ConfTrkTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

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

  ConfigurableAxis ConfV0TempFitVarBins{"ConfV0TempFitVarBins", {300, 0.95, 1.}, "V0: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfV0TempFitVarpTBins{"ConfV0TempFitVarpTBins", {20, 0.5, 4.05}, "V0: pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfV0TempFitVarInvMassBins{"ConfV0TempFitVarInvMassBins", {200, 1, 1.2}, "V0: InvMass binning"};

  ConfigurableAxis ConfV0ChildTempFitVarBins{"ConfV0ChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfV0ChildTempFitVarpTBins{"ConfV0ChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};

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

  /// Correlation part
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<bool> ConfExtendedPlots{"ConfExtendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
  Configurable<float> ConfHighkstarCut{"ConfHighkstarCut", -1., "Set a cut for high k*, above which the pairs are rejected. Set it to -1 to deactivate it"};

  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfkstarBins{"ConfkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfkTBins{"ConfkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfmTBins{"ConfmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfCPRdeltaPhiMax{"ConfCPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaMax{"ConfCPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};

  ConfigurableAxis ConfmTBins3D{"ConfmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfmultBins3D{"ConfMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfDummy{"ConfDummy", {1, 0, 1}, "Dummy axis"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

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
    trackHistoPartOne.init(&qaRegistry, ConfTrkTempFitVarpTBins, ConfTrkTempFitVarBins, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, ConfTrk1_PDGCode);
    trackHistoPartTwo.init(&qaRegistry, ConfV0TempFitVarpTBins, ConfV0TempFitVarBins, ConfDummy, ConfDummy, ConfDummy, ConfV0TempFitVarInvMassBins, ConfIsMC, ConfV02_PDGCode);
    posChildHistos.init(&qaRegistry, ConfV0ChildTempFitVarpTBins, ConfV0ChildTempFitVarBins, ConfDummy, ConfDummy, ConfDummy, ConfDummy, false, false);
    negChildHistos.init(&qaRegistry, ConfV0ChildTempFitVarpTBins, ConfV0ChildTempFitVarBins, ConfDummy, ConfDummy, ConfDummy, ConfDummy, false, false);

    sameEventCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfIsMC, ConfUse3D, ConfExtendedPlots, ConfHighkstarCut);
    sameEventCont.setPDGCodes(ConfTrk1_PDGCode, ConfV02_PDGCode);
    mixedEventCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfIsMC, ConfUse3D, ConfExtendedPlots, ConfHighkstarCut);
    mixedEventCont.setPDGCodes(ConfTrk1_PDGCode, ConfV02_PDGCode);
    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax.value, ConfCPRdeltaEtaMax.value, ConfCPRPlotPerRadii.value);
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
          MaskBit = static_cast<femtodreamcollision::BitMaskType>(mask.to_ulong());
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
  template <bool isMC, typename T, typename R, typename S>
  void doSameEvent(T& SliceTrk1, R& SliceV02, S& parts, float magFieldTesla, int multCol)
  {
    /// Histogramming same event
    for (auto& part : SliceTrk1) {
      trackHistoPartOne.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt);
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
        trackHistoPartTwo.fillQA<isMC, false>(v0, aod::femtodreamparticle::kPt);
        posChildHistos.fillQA<false, false>(posChild, aod::femtodreamparticle::kPt);
        negChildHistos.fillQA<false, false>(negChild, aod::femtodreamparticle::kPt);
      }
    }

    /// Now build the combinations
    for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceV02))) {
      const auto& posChild = parts.iteratorAt(p2.index() - 2);
      const auto& negChild = parts.iteratorAt(p2.index() - 1);
      // check cuts on V0 children
      if (((posChild.cut() & ConfV02_ChildPos_CutBit) == ConfV02_ChildPos_CutBit) &&
          ((posChild.pidcut() & ConfV02_ChildPos_TPCBit) == ConfV02_ChildPos_TPCBit) &&
          ((negChild.cut() & ConfV02_ChildNeg_CutBit) == ConfV02_ChildNeg_CutBit) &&
          ((negChild.pidcut() & ConfV02_ChildNeg_TPCBit) == ConfV02_ChildNeg_TPCBit)) {
        if (ConfIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
            continue;
          }
        }
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        sameEventCont.setPair<isMC>(p1, p2, multCol, ConfUse3D, ConfExtendedPlots);
      }
    }
  }

  void processSameEvent(o2::aod::FDCollision const& col, FilteredFDParticles const& parts)
  {
    eventHisto.fillQA(col);
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceV02 = PartitionV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceTrk1.size() == 0 && SliceV02.size() == 0) {
      return;
    }
    doSameEvent<false>(SliceTrk1, SliceV02, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEvent, "Enable processing same event", true);

  void processSameEventMasked(MaskedCollision const& col, FilteredFDParticles const& parts)
  {
    if ((col.bitmaskTrackOne() & MaskBit) != MaskBit && (col.bitmaskTrackTwo() & MaskBit) != MaskBit) {
      return;
    }
    eventHisto.fillQA(col);
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceV02 = PartitionV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<false>(SliceTrk1, SliceV02, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEventMasked, "Enable processing same event with masks", false);

  void processSameEventMC(o2::aod::FDCollision& col, FilteredFDMCParts& parts, o2::aod::FDMCParticles&)
  {
    eventHisto.fillQA(col);
    auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceMCV02 = PartitionMCV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceMCTrk1.size() == 0 && SliceMCV02.size() == 0) {
      return;
    }
    doSameEvent<true>(SliceMCTrk1, SliceMCV02, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEventMC, "Enable processing same event MC", false);

  void processSameEventMCMasked(MaskedCollision& col, FilteredFDMCParts& parts, o2::aod::FDMCParticles&)
  {
    if ((col.bitmaskTrackOne() & MaskBit) != MaskBit && (col.bitmaskTrackTwo() & MaskBit) != MaskBit) {
      return;
    }
    eventHisto.fillQA(col);
    auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceMCV02 = PartitionMCV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<true>(SliceMCTrk1, SliceMCV02, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEventMCMasked, "Enable processing same event MC with masks", false);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType SliceTrk1, PartitionType SliceV02, PartType parts, float magFieldTesla, int multCol)
  {
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
      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
      }
      // pair cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      mixedEventCont.setPair<isMC>(p1, p2, multCol, ConfUse3D, ConfExtendedPlots);
    }
  }

  void processMixedEvent(o2::aod::FDCollisions& cols, FilteredFDParticles& parts)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ConfNEventsMix, -1, cols, cols)) {
      auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceV02 = PartitionV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);

      const int multCol = collision1.multNtr();
      const auto magFieldTesla1 = collision1.magField();
      if (SliceTrk1.size() == 0 || SliceV02.size() == 0) {
        return;
      }
      doMixedEvent<false>(SliceTrk1, SliceV02, parts, magFieldTesla1, multCol);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processMixedEvent, "Enable processing mixed events", true);

  void processMixedEventMasked(MaskedCollisions const& cols, FilteredFDParticles const& parts)
  {
    Partition<MaskedCollisions> PartitionMaskedCol1 = (aod::femtodreamcollision::bitmaskTrackOne & MaskBit) == MaskBit;
    Partition<MaskedCollisions> PartitionMaskedCol2 = (aod::femtodreamcollision::bitmaskTrackTwo & MaskBit) == MaskBit;

    PartitionMaskedCol1.bindTable(cols);
    PartitionMaskedCol2.bindTable(cols);

    for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockUpperIndexPolicy(colBinning, ConfNEventsMix.value, -1, PartitionMaskedCol1, PartitionMaskedCol2))) {
      auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceV02 = PartitionV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);

      const int multCol = collision1.multNtr();
      const auto magFieldTesla1 = collision1.magField();
      doMixedEvent<false>(SliceTrk1, SliceV02, parts, magFieldTesla1, multCol);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processMixedEventMasked, "Enable processing mixed events with masks", false);

  void processMixedEventMC(o2::aod::FDCollisions& cols, FilteredFDMCParts& parts, o2::aod::FDMCParticles&)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ConfNEventsMix, -1, cols, cols)) {
      auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceMCV02 = PartitionMCV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);

      const int multCol = collision1.multNtr();
      const auto magFieldTesla1 = collision1.magField();
      if (SliceMCTrk1.size() == 0 || SliceMCV02.size() == 0) {
        return;
      }
      doMixedEvent<true>(SliceMCTrk1, SliceMCV02, parts, magFieldTesla1, multCol);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processMixedEventMC, "Enable processing mixed events MC", false);

  void processMixedEventMCMasked(MaskedCollisions& cols, FilteredFDMCParts& parts, o2::aod::FDMCParticles&)
  {
    Partition<MaskedCollisions> PartitionMaskedCol1 = (aod::femtodreamcollision::bitmaskTrackOne & MaskBit) == MaskBit;
    Partition<MaskedCollisions> PartitionMaskedCol2 = (aod::femtodreamcollision::bitmaskTrackTwo & MaskBit) == MaskBit;

    PartitionMaskedCol1.bindTable(cols);
    PartitionMaskedCol2.bindTable(cols);

    for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockUpperIndexPolicy(colBinning, ConfNEventsMix.value, -1, PartitionMaskedCol1, PartitionMaskedCol2))) {
      auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceMCV02 = PartitionMCV02->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      const int multCol = collision1.multNtr();
      const auto magFieldTesla1 = collision1.magField();
      doMixedEvent<true>(SliceMCTrk1, SliceMCV02, parts, magFieldTesla1, multCol);
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
