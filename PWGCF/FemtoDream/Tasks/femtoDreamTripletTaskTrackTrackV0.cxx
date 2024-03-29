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

/// \file femtoDreamTripletTaskTrackTrackV0.cxx
/// \brief Tasks that reads the track and V0 tables and creates triplets; only two identical tracks and a V0 can be used
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de

#include <vector>
#include <bitset>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/Core/femtoDreamContainerThreeBody.h"
#include "PWGCF/FemtoDream/Core/femtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoDreamTripletTaskTrackTrackV0 {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  using MaskedCollisions = soa::Join<aod::FDCollisions, aod::FDColMasks>;
  using MaskedCollision = MaskedCollisions::iterator;
  aod::femtodreamcollision::BitMaskType MaskBit = -1;

  Configurable<bool> ConfMixIfTripletPresent{"ConfMixIfTripletPresent", true, "Use for mixing only events which have a TTV0 triplet"};
  Configurable<bool> ConfMixIfTVOPairPresent{"ConfMixIfTVOPairPresent", false, "Use for mixing only events which have a TV0 pair (at least one track and one V0)"};
  Configurable<bool> ConfMixIfTOrVOPartsPresent{"ConfMixIfTOrVOPartsPresent", false, "Use for mixing only events which have at least one particle of interest"};

  /// Track selection
  Configurable<int> ConfPDGCodePart{"ConfPDGCodePart", 2212, "Particle PDG code"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfCutPart{"ConfCutPart", 5542474, "Track - Selection bit from cutCulator"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfTPCPIDBit{"ConfTPCPIDBit", 16, "PID TPC bit from cutCulator "};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfTPCTOFPIDBit{"ConfTPCTOFPIDBit", 8, "PID TPCTOF bit from cutCulator"};
  Configurable<float> ConfMaxpT{"ConfMaxpT", 4.05f, "Maximum transverse momentum of the particles"};
  Configurable<float> ConfPIDthrMom{"ConfPIDthrMom", 1.f, "Momentum threshold from which TPC and TOF are required for PID"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};

  /// Partition for selected particles
  Partition<aod::FDParticles> SelectedParts = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                              ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfPIDthrMom, ncheckbit(aod::femtodreamparticle::pidcut, ConfTPCPIDBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTPCTOFPIDBit)) &&
                                              (ncheckbit(aod::femtodreamparticle::cut, ConfCutPart)) &&
                                              (aod::femtodreamparticle::pt < ConfMaxpT);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> SelectedPartsMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                                            ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfPIDthrMom, ncheckbit(aod::femtodreamparticle::pidcut, ConfTPCPIDBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTPCTOFPIDBit)) &&
                                                                            (ncheckbit(aod::femtodreamparticle::cut, ConfCutPart)) &&
                                                                            (aod::femtodreamparticle::pt < ConfMaxpT);

  /// Histogramming of selected tracks
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoSelectedParts;

  /// V0 selection
  Configurable<int> ConfPDGCodeV0{"ConfPDGCodeV0", 3122, "V0 PDG code"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfCutV0{"ConfCutV0", 338, "V0 - Selection bit from cutCulator"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> Conf_ChildPos_CutV0{"Conf_ChildPos_CutV0", 149, "Selection bit for positive child of V0"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> Conf_ChildPos_TPCBitV0{"Conf_ChildPos_TPCBitV0", 2, "PID TPC bit for positive child of V0"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> Conf_ChildNeg_CutV0{"Conf_ChildNeg_CutV0", 149, "Selection bit for negative child of V0"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> Conf_ChildNeg_TPCBitV0{"Conf_ChildNeg_TPCBitV0", 2, "PID TPC bit for negative child of V0"};

  Configurable<float> Conf_minInvMass_V0{"Conf_minInvMass_V0", 1.08, "Minimum invariant mass of V0 (particle)"};
  Configurable<float> Conf_maxInvMass_V0{"Conf_maxInvMass_V0", 1.15, "Maximum invariant mass of V0 (particle)"};
  Configurable<float> Conf_minInvMassAnti_V0{"Conf_minInvMassAnti_V0", 1.08, "Minimum invariant mass of V0 (antiparticle)"};
  Configurable<float> Conf_maxInvMassAnti_V0{"Conf_maxInvMassAnti_V0", 1.15, "Maximum invariant mass of V0 (antiparticle)"};

  Configurable<float> Conf_minPt_V0{"Conf_minPt_V0", 0., "Minimum pT of V0"};
  Configurable<float> Conf_maxPt_V0{"Conf_maxPt_V0", 999., "Maximum pT of V0"};

  // Partition for selected particles
  Partition<aod::FDParticles> SelectedV0s = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                            (ncheckbit(aod::femtodreamparticle::cut, ConfCutV0)) &&
                                            (aod::femtodreamparticle::mLambda > Conf_minInvMass_V0) &&
                                            (aod::femtodreamparticle::mLambda < Conf_maxInvMass_V0) &&
                                            (aod::femtodreamparticle::mAntiLambda > Conf_minInvMassAnti_V0) &&
                                            (aod::femtodreamparticle::mAntiLambda < Conf_maxInvMassAnti_V0) &&
                                            (aod::femtodreamparticle::pt > Conf_minPt_V0) &&
                                            (aod::femtodreamparticle::pt < Conf_maxPt_V0);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> SelectedV0sMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                                                          (ncheckbit(aod::femtodreamparticle::cut, ConfCutV0)) &&
                                                                          (aod::femtodreamparticle::mLambda > Conf_minInvMass_V0) &&
                                                                          (aod::femtodreamparticle::mLambda < Conf_maxInvMass_V0) &&
                                                                          (aod::femtodreamparticle::mAntiLambda > Conf_minInvMassAnti_V0) &&
                                                                          (aod::femtodreamparticle::mAntiLambda < Conf_maxInvMassAnti_V0) &&
                                                                          (aod::femtodreamparticle::pt > Conf_minPt_V0) &&
                                                                          (aod::femtodreamparticle::pt < Conf_maxPt_V0);

  /// Histogramming of selected V0s
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 2> particleHistoSelectedV0s;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 3> particleHistoPosChild;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 4> particleHistoNegChild;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// Particles
  ConfigurableAxis ConfTempFitVarBinsTrack{"ConfTempFitVarBinsTrack", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (track)"};
  ConfigurableAxis ConfTempFitVarBinsV0{"ConfTempFitVarBinsV0", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis ConfTempFitVarBinsV0Child{"ConfTempFitVarBinsV0Child", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0 child)"};
  ConfigurableAxis ConfTempFitVarpTV0Bins{"ConfTempFitVarpTV0Bins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis ConfTempFitVarpTV0Child{"ConfTempFitVarpTV0Child", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0 child)"};
  ConfigurableAxis ConfInvMassBins{"ConfInvMassBins", {200, 1, 1.2}, "InvMass binning"};

  /// Correlations
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

  ConfigurableAxis ConfQ3Bins{"ConfQ3Bins", {2000, 0., 8.}, "binning Q3"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfCPRdeltaPhiMax{"ConfCPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaMax{"ConfCPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
  ConfigurableAxis ConfDummy{"ConfDummy", {1, 0, 1}, "Dummy axis"};

  FemtoDreamContainerThreeBody<femtoDreamContainerThreeBody::EventType::same, femtoDreamContainerThreeBody::Observable::Q3> sameEventCont;
  FemtoDreamContainerThreeBody<femtoDreamContainerThreeBody::EventType::mixed, femtoDreamContainerThreeBody::Observable::Q3> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleanerTrackTrack;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCleanerTrackV0;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejectionTrackTrack;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCloseRejectionTrackV0;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry ThreeBodyQARegistry{"ThreeBodyQARegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext& context)
  {

    eventHisto.init(&qaRegistry);
    trackHistoSelectedParts.init(&qaRegistry, ConfDummy, ConfDummy, ConfTempFitVarpTBins, ConfDummy, ConfDummy, ConfTempFitVarBinsTrack, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, ConfPDGCodePart);
    particleHistoSelectedV0s.init(&qaRegistry, ConfDummy, ConfDummy, ConfTempFitVarpTV0Bins, ConfDummy, ConfDummy, ConfTempFitVarBinsV0, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfInvMassBins, ConfIsMC, ConfPDGCodeV0);
    particleHistoPosChild.init(&qaRegistry, ConfDummy, ConfDummy, ConfTempFitVarpTV0Child, ConfDummy, ConfDummy, ConfTempFitVarBinsV0Child, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, 0);
    particleHistoNegChild.init(&qaRegistry, ConfDummy, ConfDummy, ConfTempFitVarpTV0Child, ConfDummy, ConfDummy, ConfTempFitVarBinsV0Child, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, 0);

    ThreeBodyQARegistry.add("TripletTaskQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    ThreeBodyQARegistry.add("TripletTaskQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    std::vector<double> tmpVecMult = ConfMultBins;
    framework::AxisSpec multAxis = {tmpVecMult, "Multiplicity"};
    ThreeBodyQARegistry.add("TripletTaskQA/hSEMultVSGoodTracks", ";Mult;GoodT", kTH2F, {multAxis, {100, 0, 100}});

    sameEventCont.init(&resultRegistry, ConfQ3Bins, ConfMultBins, ConfIsMC);
    mixedEventCont.init(&resultRegistry, ConfQ3Bins, ConfMultBins, ConfIsMC);
    sameEventCont.setPDGCodes(ConfPDGCodePart, ConfPDGCodePart, ConfPDGCodeV0);
    mixedEventCont.setPDGCodes(ConfPDGCodePart, ConfPDGCodePart, ConfPDGCodeV0);

    pairCleanerTrackTrack.init(&qaRegistry);
    pairCleanerTrackV0.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejectionTrackTrack.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax.value, ConfCPRdeltaEtaMax.value, ConfCPRPlotPerRadii.value);
      pairCloseRejectionTrackV0.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax.value, ConfCPRdeltaEtaMax.value, ConfCPRPlotPerRadii.value);
    }

    // get bit for the collision mask
    std::bitset<8 * sizeof(aod::femtodreamcollision::BitMaskType)> mask;
    int index = 0;
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.find("femto-dream-triplet-task-track-track-v0") != std::string::npos) {
        if (containsNameValuePair(device.options, "ConfCutPart", ConfCutPart.value) &&
            containsNameValuePair(device.options, "ConfTPCPIDBit", ConfTPCPIDBit.value) &&
            containsNameValuePair(device.options, "ConfTPCTOFPIDBit", ConfTPCTOFPIDBit.value) &&
            containsNameValuePair(device.options, "ConfPIDthrMom", ConfPIDthrMom.value) &&
            containsNameValuePair(device.options, "ConfMaxpT", ConfMaxpT.value) &&
            containsNameValuePair(device.options, "ConfCutV0", ConfCutV0.value) &&
            containsNameValuePair(device.options, "Conf_ChildPos_CutV0", Conf_ChildPos_CutV0.value) &&
            containsNameValuePair(device.options, "Conf_ChildPos_TPCBitV0", Conf_ChildPos_TPCBitV0.value) &&
            containsNameValuePair(device.options, "Conf_ChildNeg_CutV0", Conf_ChildNeg_CutV0.value) &&
            containsNameValuePair(device.options, "Conf_ChildNeg_TPCBitV0", Conf_ChildNeg_TPCBitV0.value) &&
            containsNameValuePair(device.options, "Conf_minInvMass_V0", Conf_minInvMass_V0.value) &&
            containsNameValuePair(device.options, "Conf_maxInvMass_V0", Conf_maxInvMass_V0.value) &&
            containsNameValuePair(device.options, "Conf_minInvMassAnti_V0", Conf_minInvMassAnti_V0.value) &&
            containsNameValuePair(device.options, "Conf_maxInvMassAnti_V0", Conf_maxInvMassAnti_V0.value) &&
            containsNameValuePair(device.options, "Conf_minPt_V0", Conf_minPt_V0.value) &&
            containsNameValuePair(device.options, "Conf_maxPt_V0", Conf_maxPt_V0.value)) {
          mask.set(index);
          MaskBit = static_cast<aod::femtodreamcollision::BitMaskType>(mask.to_ulong());
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

    if ((ConfMixIfTripletPresent && ConfMixIfTVOPairPresent) ||
        (ConfMixIfTripletPresent && ConfMixIfTOrVOPartsPresent) ||
        (ConfMixIfTVOPairPresent && ConfMixIfTOrVOPartsPresent)) {
      LOG(fatal) << "Only one method of mixing can be chosen!";
    }
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
    eventHisto.fillQA(col);
  }

  /// This function processes the same event and takes care of all the histogramming
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupSelectedTracks partition for the first particle passed by the process function
  /// @param parts femtoDreamParticles table (in case of Monte Carlo joined with FemtoDreamMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doSameEvent(PartitionType groupSelectedTracks, PartitionType groupSelectedV0s, PartType parts, float magFieldTesla, int multCol, float centCol)
  {
    /// Histograming tracks
    for (auto& part : groupSelectedTracks) {
      trackHistoSelectedParts.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt, multCol, centCol);
    }
    /// Histograming V0s
    for (auto& V0 : groupSelectedV0s) {
      const auto& posChild = parts.iteratorAt(V0.index() - 2);
      const auto& negChild = parts.iteratorAt(V0.index() - 1);

      particleHistoSelectedV0s.fillQA<isMC, false>(V0, aod::femtodreamparticle::kPt, multCol, centCol);
      particleHistoPosChild.fillQA<isMC, false>(posChild, aod::femtodreamparticle::kPt, multCol, centCol);
      particleHistoNegChild.fillQA<isMC, false>(negChild, aod::femtodreamparticle::kPt, multCol, centCol);

      if (((posChild.cut() & Conf_ChildPos_CutV0) == Conf_ChildPos_CutV0 &&
           (posChild.pidcut() & Conf_ChildPos_TPCBitV0) == Conf_ChildPos_TPCBitV0 &&
           (negChild.cut() & Conf_ChildNeg_CutV0) == Conf_ChildNeg_CutV0 &&
           (negChild.pidcut() & Conf_ChildNeg_TPCBitV0) == Conf_ChildNeg_TPCBitV0)) {
        particleHistoSelectedV0s.fillQA<isMC, false>(V0, aod::femtodreamparticle::kPt, multCol, centCol);
        particleHistoPosChild.fillQA<isMC, false>(posChild, aod::femtodreamparticle::kPt, multCol, centCol);
        particleHistoNegChild.fillQA<isMC, false>(negChild, aod::femtodreamparticle::kPt, multCol, centCol);
      }
    }

    /// Now build the combinations
    for (auto& V0 : groupSelectedV0s) {
      const auto& posChild = parts.iteratorAt(V0.index() - 2);
      const auto& negChild = parts.iteratorAt(V0.index() - 1);
      if (!((posChild.cut() & Conf_ChildPos_CutV0) == Conf_ChildPos_CutV0 &&
            (posChild.pidcut() & Conf_ChildPos_TPCBitV0) == Conf_ChildPos_TPCBitV0 &&
            (negChild.cut() & Conf_ChildNeg_CutV0) == Conf_ChildNeg_CutV0 &&
            (negChild.pidcut() & Conf_ChildNeg_TPCBitV0) == Conf_ChildNeg_TPCBitV0)) {
        continue;
      }

      for (auto& [T1, T2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupSelectedTracks, groupSelectedTracks))) {

        // Close pair rejection
        if (ConfIsCPR.value) {
          if (pairCloseRejectionTrackTrack.isClosePair(T1, T2, parts, magFieldTesla)) {
            continue;
          }
          if (pairCloseRejectionTrackV0.isClosePair(T1, V0, parts, magFieldTesla)) {
            continue;
          }
          if (pairCloseRejectionTrackV0.isClosePair(T2, V0, parts, magFieldTesla)) {
            continue;
          }
        }

        // track cleaning
        if (!pairCleanerTrackTrack.isCleanPair(T1, T2, parts)) {
          continue;
        }
        if (!pairCleanerTrackV0.isCleanPair(T2, V0, parts)) {
          continue;
        }
        if (!pairCleanerTrackV0.isCleanPair(T1, V0, parts)) {
          continue;
        }
        sameEventCont.setTriplet<isMC>(T1, T2, V0, multCol);
      }
    }
  }

  /// process function to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEvent(o2::aod::FDCollision& col,
                        o2::aod::FDParticles& parts)
  {
    fillCollision(col);
    auto thegroupSelectedTracks = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedV0s = SelectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (thegroupSelectedTracks.size() < 2 || thegroupSelectedV0s.size() < 1) {
      return;
    }
    doSameEvent<false>(thegroupSelectedTracks, thegroupSelectedV0s, parts, col.magField(), col.multNtr(), col.multV0M());
  }
  PROCESS_SWITCH(femtoDreamTripletTaskTrackTrackV0, processSameEvent, "Enable processing same event", true);

  /// process function to call doSameEvent with Data which has a mask for containing particles or not
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEventMasked(MaskedCollision& col, o2::aod::FDParticles& parts)
  {
    fillCollision(col);
    auto thegroupSelectedTracks = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedV0s = SelectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (thegroupSelectedTracks.size() < 2 || thegroupSelectedV0s.size() < 1) {
      return;
    }
    doSameEvent<false>(thegroupSelectedTracks, thegroupSelectedV0s, parts, col.magField(), col.multNtr(), col.multV0M());
  }
  PROCESS_SWITCH(femtoDreamTripletTaskTrackTrackV0, processSameEventMasked, "Enable processing same event with masks", false);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(o2::aod::FDCollision& col,
                          soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                          o2::aod::FDMCParticles&)
  {
    fillCollision(col);
    auto thegroupSelectedTracks = SelectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedV0s = SelectedV0sMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (thegroupSelectedTracks.size() < 2 || thegroupSelectedV0s.size() < 1) {
      return;
    }
    doSameEvent<true>(thegroupSelectedTracks, thegroupSelectedV0s, parts, col.magField(), col.multNtr(), col.multV0M());
  }
  PROCESS_SWITCH(femtoDreamTripletTaskTrackTrackV0, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  /// process function for to call doSameEvent with Monte Carlo which has a mask for containing particles or not
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMCMasked(MaskedCollision& col,
                                soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                                o2::aod::FDMCParticles&)
  {
    fillCollision(col);
    auto thegroupSelectedTracks = SelectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedV0s = SelectedV0sMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (thegroupSelectedTracks.size() < 2 || thegroupSelectedV0s.size() < 1) {
      return;
    }
    doSameEvent<true>(thegroupSelectedTracks, thegroupSelectedV0s, parts, col.magField(), col.multNtr(), col.multV0M());
  }
  PROCESS_SWITCH(femtoDreamTripletTaskTrackTrackV0, processSameEventMCMasked, "Enable processing same event for Monte Carlo", false);

  /// This function processes the mixed event
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsOne partition for the first particle passed by the process function
  /// \param groupPartsTwo partition for the second particle passed by the process function
  /// \param groupPartsThree partition for the third particle passed by the process function
  /// \param parts femtoDreamParticles table (in case of Monte Carlo joined with FemtoDreamMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartitionType groupPartsThree, PartType parts, float magFieldTesla, int multCol)
  {
    for (auto& [T1, T2, V0] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo, groupPartsThree))) {
      const auto& posChild = parts.iteratorAt(V0.globalIndex() - 2);
      const auto& negChild = parts.iteratorAt(V0.globalIndex() - 1);

      if (!((posChild.cut() & Conf_ChildPos_CutV0) == Conf_ChildPos_CutV0 &&
            (posChild.pidcut() & Conf_ChildPos_TPCBitV0) == Conf_ChildPos_TPCBitV0 &&
            (negChild.cut() & Conf_ChildNeg_CutV0) == Conf_ChildNeg_CutV0 &&
            (negChild.pidcut() & Conf_ChildNeg_TPCBitV0) == Conf_ChildNeg_TPCBitV0)) {
        continue;
      }

      // Close pair rejection
      if (ConfIsCPR.value) {
        if (pairCloseRejectionTrackTrack.isClosePair(T1, T2, parts, magFieldTesla)) {
          continue;
        }
        if (pairCloseRejectionTrackV0.isClosePair(T1, V0, parts, magFieldTesla)) {
          continue;
        }
        if (pairCloseRejectionTrackV0.isClosePair(T2, V0, parts, magFieldTesla)) {
          continue;
        }
      }
      mixedEventCont.setTriplet<isMC>(T1, T2, V0, multCol);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoDreamParticleTable
  void processMixedEvent(o2::aod::FDCollisions& cols,
                         o2::aod::FDParticles& parts)
  {
    for (auto& [collision1, collision2, collision3] : soa::selfCombinations(colBinning, ConfNEventsMix, -1, cols, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsThree = SelectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();
      const auto& magFieldTesla3 = collision3.magField();

      if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
        continue;
      }
      // CONSIDER testing different strategies to which events to use

      doMixedEvent<false>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoDreamTripletTaskTrackTrackV0, processMixedEvent, "Enable processing mixed events", true);

  /// process function for to call doMixedEvent with Data which has a mask for containing particles or not
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoDreamParticleTable
  void processMixedEventMasked(MaskedCollisions& cols, o2::aod::FDParticles& parts)
  {
    if (ConfMixIfTripletPresent || ConfMixIfTVOPairPresent) {
      Partition<MaskedCollisions> PartitionMaskedCol1 = (ConfMixIfTripletPresent && (aod::femtodreamcollision::bitmaskTrackTwo & MaskBit) == MaskBit && (aod::femtodreamcollision::bitmaskTrackThree & MaskBit) == MaskBit) ||
                                                        (ConfMixIfTVOPairPresent && (aod::femtodreamcollision::bitmaskTrackOne & MaskBit) == MaskBit && (aod::femtodreamcollision::bitmaskTrackThree & MaskBit) == MaskBit);
      PartitionMaskedCol1.bindTable(cols);

      for (auto& [collision1, collision2, collision3] : soa::selfCombinations(colBinning, ConfNEventsMix, -1, *PartitionMaskedCol1.mFiltered, *PartitionMaskedCol1.mFiltered, *PartitionMaskedCol1.mFiltered)) {

        const int multiplicityCol = collision1.multNtr();
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

        auto groupPartsOne = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
        auto groupPartsThree = SelectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

        const auto& magFieldTesla1 = collision1.magField();
        const auto& magFieldTesla2 = collision2.magField();
        const auto& magFieldTesla3 = collision3.magField();

        if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
          continue;
        }
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
      }
    } else if (ConfMixIfTOrVOPartsPresent) {
      Partition<MaskedCollisions> PartitionMaskedColT = ((aod::femtodreamcollision::bitmaskTrackOne & MaskBit) == MaskBit);
      Partition<MaskedCollisions> PartitionMaskedColV0 = ((aod::femtodreamcollision::bitmaskTrackThree & MaskBit) == MaskBit);
      PartitionMaskedColT.bindTable(cols);
      PartitionMaskedColV0.bindTable(cols);

      for (auto& [ColWithTrack1, ColWithTrack2, ColWithV0] : soa::combinations(soa::CombinationsBlockFullIndexPolicy(colBinning, ConfNEventsMix, -1, *PartitionMaskedColT.mFiltered, *PartitionMaskedColT.mFiltered, *PartitionMaskedColV0.mFiltered))) {
        if (ColWithTrack1.globalIndex() == ColWithTrack2.globalIndex() || ColWithTrack1.globalIndex() == ColWithV0.globalIndex() || ColWithTrack2.globalIndex() == ColWithV0.globalIndex()) {
          continue;
        }
        if (ColWithTrack1.globalIndex() > ColWithTrack2.globalIndex()) {
          continue;
        }
        const int multiplicityCol = ColWithTrack1.multNtr();
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({ColWithTrack1.posZ(), multiplicityCol}));

        auto groupPartsOne = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithTrack1.globalIndex(), cache);
        auto groupPartsTwo = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithTrack2.globalIndex(), cache);
        auto groupPartsThree = SelectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithV0.globalIndex(), cache);

        const auto& magFieldTesla1 = ColWithTrack1.magField();
        const auto& magFieldTesla2 = ColWithTrack2.magField();
        const auto& magFieldTesla3 = ColWithV0.magField();

        if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
          continue;
        }
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
      }
    }
  }
  PROCESS_SWITCH(femtoDreamTripletTaskTrackTrackV0, processMixedEventMasked, "Enable processing mixed events", false);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(o2::aod::FDCollisions& cols,
                           soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                           o2::aod::FDMCParticles&)
  {
    for (auto& [collision1, collision2, collision3] : soa::selfCombinations(colBinning, ConfNEventsMix, -1, cols, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = SelectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = SelectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsThree = SelectedV0sMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();
      const auto& magFieldTesla3 = collision3.magField();

      if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
        continue;
      }
      // CONSIDER testing different strategies to which events to use

      doMixedEvent<true>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoDreamTripletTaskTrackTrackV0, processMixedEventMC, "Enable processing mixed events MC", false);

  /// brief process function for to call doMixedEvent with Monte Carlo which has a mask for containing particles or not
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMCMasked(MaskedCollisions& cols,
                                 soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                                 o2::aod::FDMCParticles&)
  {
    if (ConfMixIfTripletPresent || ConfMixIfTVOPairPresent) {
      Partition<MaskedCollisions> PartitionMaskedCol1 = (ConfMixIfTripletPresent && (aod::femtodreamcollision::bitmaskTrackTwo & MaskBit) == MaskBit && (aod::femtodreamcollision::bitmaskTrackThree & MaskBit) == MaskBit) ||
                                                        (ConfMixIfTVOPairPresent && (aod::femtodreamcollision::bitmaskTrackOne & MaskBit) == MaskBit && (aod::femtodreamcollision::bitmaskTrackThree & MaskBit) == MaskBit);
      PartitionMaskedCol1.bindTable(cols);

      for (auto& [collision1, collision2, collision3] : soa::selfCombinations(colBinning, ConfNEventsMix, -1, *PartitionMaskedCol1.mFiltered, *PartitionMaskedCol1.mFiltered, *PartitionMaskedCol1.mFiltered)) {
        const int multiplicityCol = collision1.multNtr();
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

        auto groupPartsOne = SelectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = SelectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
        auto groupPartsThree = SelectedV0sMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

        const auto& magFieldTesla1 = collision1.magField();
        const auto& magFieldTesla2 = collision2.magField();
        const auto& magFieldTesla3 = collision3.magField();

        if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
          continue;
        }
        doMixedEvent<true>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
      }
    } else if (ConfMixIfTOrVOPartsPresent) {
      Partition<MaskedCollisions> PartitionMaskedColT = ((aod::femtodreamcollision::bitmaskTrackOne & MaskBit) == MaskBit);
      Partition<MaskedCollisions> PartitionMaskedColV0 = ((aod::femtodreamcollision::bitmaskTrackThree & MaskBit) == MaskBit);
      PartitionMaskedColT.bindTable(cols);
      PartitionMaskedColV0.bindTable(cols);

      for (auto& [ColWithTrack1, ColWithTrack2, ColWithV0] : soa::combinations(soa::CombinationsBlockFullIndexPolicy(colBinning, ConfNEventsMix, -1, *PartitionMaskedColT.mFiltered, *PartitionMaskedColT.mFiltered, *PartitionMaskedColV0.mFiltered))) {
        if (ColWithTrack1.globalIndex() == ColWithTrack2.globalIndex() || ColWithTrack1.globalIndex() == ColWithV0.globalIndex() || ColWithTrack2.globalIndex() == ColWithV0.globalIndex()) {
          continue;
        }
        if (ColWithTrack1.globalIndex() > ColWithTrack2.globalIndex()) {
          continue;
        }
        const int multiplicityCol = ColWithTrack1.multNtr();
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({ColWithTrack1.posZ(), multiplicityCol}));

        auto groupPartsOne = SelectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithTrack1.globalIndex(), cache);
        auto groupPartsTwo = SelectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithTrack2.globalIndex(), cache);
        auto groupPartsThree = SelectedV0sMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithV0.globalIndex(), cache);

        const auto& magFieldTesla1 = ColWithTrack1.magField();
        const auto& magFieldTesla2 = ColWithTrack2.magField();
        const auto& magFieldTesla3 = ColWithV0.magField();

        if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
          continue;
        }
        doMixedEvent<true>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
      }
    }
  }
  PROCESS_SWITCH(femtoDreamTripletTaskTrackTrackV0, processMixedEventMCMasked, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamTripletTaskTrackTrackV0>(cfgc),
  };
  return workflow;
}
