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

/// \file femtoDreamTripletTaskTrackTrackCascade.cxx
/// \brief Tasks that reads the track and Cascade tables and creates triplets; only two identical tracks and a Cascade can be used
/// \author Raffaele Del Grande, CTU Prague

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamContainerThreeBody.h"
#include "PWGCF/FemtoDream/Core/femtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include "TDatabasePDG.h"

#include <bitset>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoDreamTripletTaskTrackTrackCascade {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  using MaskedCollisions = soa::Join<aod::FDCollisions, aod::FDColMasks>;
  using MaskedCollision = MaskedCollisions::iterator;
  aod::femtodreamcollision::BitMaskType MaskBit = -1;
  float mMassOne = -999, mMassTwo = -999, mMassThree = -999;

  Configurable<bool> ConfMixIfTripletPresent{"ConfMixIfTripletPresent", true, "Use for mixing only events which have a TTV0 triplet"};
  Configurable<bool> ConfMixIfTVOPairPresent{"ConfMixIfTVOPairPresent", false, "Use for mixing only events which have a TV0 pair (at least one track and one V0)"};
  Configurable<bool> ConfMixIfTOrVOPartsPresent{"ConfMixIfTOrVOPartsPresent", false, "Use for mixing only events which have at least one particle of interest"};

  // which CPR to use, old is with a possible bug and new is fixed
  Configurable<bool> ConfUseOLD_possiblyWrong_CPR{"ConfUseOLD_possiblyWrong_CPR", true, "Use for old CPR, which possibly has a bug. This is implemented only for debugging reasons to compare old and new code on hyperloop datasets."};
  Configurable<int> ConfAtWhichRadiiToCut{"ConfAtWhichRadiiToCut", 1, "At which radii perform deta dphi selection: 0 - at PV, 1 - averaged phi, 2 - at given radii"};
  Configurable<float> ConfAtWhichTPCRadii{"ConfAtWhichTPCRadii", 85., "If ConfAtWhichRadiiToCut = 2; this allows to select at which TPC radii to cut"};

  /// Track selection
  Configurable<int> ConfPDGCodePart{"ConfPDGCodePart", 2212, "Particle PDG code"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfCutPart{"ConfCutPart", 5542474, "Track - Selection bit from cutCulator"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfTPCPIDBit{"ConfTPCPIDBit", 16, "PID TPC bit from cutCulator "};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfTPCTOFPIDBit{"ConfTPCTOFPIDBit", 8, "PID TPCTOF bit from cutCulator"};
  Configurable<float> ConfMaxpT{"ConfMaxpT", 4.05f, "Maximum transverse momentum of the particles"};
  Configurable<float> ConfMinpT{"ConfMinpT", 0.3f, "Minimum transverse momentum of the particles"};
  Configurable<float> ConfMaxDCAxy{"ConfMaxDCAxy", -0.1f, "Maximum DCAxy of the particles"};
  Configurable<float> ConfMinDCAxy{"ConfMinDCAxy", 0.1f, "Minimum DCAxy of the particles"};
  Configurable<float> ConfPIDthrMom{"ConfPIDthrMom", 1.f, "Momentum threshold from which TPC and TOF are required for PID"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<bool> ConfDCACutPtDep{"ConfDCACutPtDep", false, "Use pt dependent dca cut for tracks"};

  /// Partition for selected particles
  Partition<aod::FDParticles> SelectedParts = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                              ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfPIDthrMom, ncheckbit(aod::femtodreamparticle::pidcut, ConfTPCPIDBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTPCTOFPIDBit)) &&
                                              (ncheckbit(aod::femtodreamparticle::cut, ConfCutPart)) &&
                                              (aod::femtodreamparticle::pt < ConfMaxpT) &&
                                              (aod::femtodreamparticle::pt > ConfMinpT) &&
                                              ifnode(ConfDCACutPtDep, (nabs(aod::femtodreamparticle::tempFitVar) <= 0.0105f + (0.035f / npow(aod::femtodreamparticle::pt, 1.1f))),
                                                     ((aod::femtodreamparticle::tempFitVar >= ConfMinDCAxy) &&
                                                      (aod::femtodreamparticle::tempFitVar <= ConfMaxDCAxy)));
  ;

  /// Histogramming of selected tracks
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoSelectedParts;

  /// Cascade selection
  Configurable<int> ConfPDGCodeCascade{"ConfPDGCodeCascade", 3312, "Cascade PDG code"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfCutCascade{"ConfCutCascade", 32221874, "Selection bit from cutCulator"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ChildPos_CutBit{"ChildPos_CutBit", 278, "Selection bit for positive child of Cascade"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ChildPos_TPCBit{"ChildPos_TPCBit", 1024, "PID TPC bit for positive child of Cascade"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ChildNeg_CutBit{"ChildNeg_CutBit", 277, "Selection bit for negative child of Cascade"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ChildNeg_TPCBit{"ChildNeg_TPCBit", 4096, "PID TPC bit for negative child of Cascade"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ChildBach_CutBit{"ChildBach_CutBit", 277, "Selection bit for negative child of Cascade"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ChildBach_TPCBit{"ChildBach_TPCBit", 64, "PID TPC bit for bachelor child of Cascade"}; // RAFFA: Why this pion has a different TPC bit? // BACHELOR means pion from the Xi

  Configurable<float> Conf_minInvMass_Cascade{"Conf_minInvMass_Cascade", 1.2, "Minimum invariant mass of Cascade (particle)"};
  Configurable<float> Conf_maxInvMass_Cascade{"Conf_maxInvMass_Cascade", 1.4, "Maximum invariant mass of Cascade (particle)"};
  Configurable<float> InvMassV0DaughMin{"InvMassV0DaughMin", 0., "Minimum invariant mass of the V0 Daughter"};
  Configurable<float> InvMassV0DaughMax{"InvMassV0DaughMax", 999., "Maximum invariant mass of the V0 Daughter"};

  Configurable<float> Conf_minPt_Cascade{"Conf_minPt_Cascade", 0., "Minimum pT of Cascade"};
  Configurable<float> Conf_maxPt_Cascade{"Conf_maxPt_Cascade", 999., "Maximum pT of Cascade"};

  // Partition for selected particles
  Partition<aod::FDParticles> SelectedCascades = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kCascade)) &&
                                                 (ncheckbit(aod::femtodreamparticle::cut, ConfCutCascade)) &&
                                                 (aod::femtodreamparticle::mLambda > Conf_minInvMass_Cascade) &&
                                                 (aod::femtodreamparticle::mLambda < Conf_maxInvMass_Cascade) &&
                                                 (aod::femtodreamparticle::mAntiLambda > InvMassV0DaughMin) &&
                                                 (aod::femtodreamparticle::mAntiLambda < InvMassV0DaughMax) &&
                                                 (aod::femtodreamparticle::pt > Conf_minPt_Cascade) &&
                                                 (aod::femtodreamparticle::pt < Conf_maxPt_Cascade);

  /// Histogramming of selected Cascades
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascade, 2> particleHistoSelectedCascades;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeV0Child, 3> particleHistoPosChild;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeV0Child, 4> particleHistoNegChild;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeBachelor, 8> particleHistoBachChild;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// Binning
  ConfigurableAxis ConfTempFitVarBinsTrack{"ConfTempFitVarBinsTrack", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (track)"};
  ConfigurableAxis ConfTempFitVarBinsCascade{"ConfTempFitVarBinsCascade", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis ConfTempFitVarBinsCascadeChild{"ConfTempFitVarBinsCascadeChild", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0 child)"};
  ConfigurableAxis ConfTempFitVarpTV0Bins{"ConfTempFitVarpTV0Bins", {35, 0, 6}, "pT binning of the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis ConfTempFitVarpTV0Child{"ConfTempFitVarpTV0Child", {35, 0, 6}, "pT binning of the pT vs. TempFitVar plot (V0 child)"};
  ConfigurableAxis ConfInvMassBins{"ConfInvMassBins", {200, 1, 1.6}, "InvMass binning"};
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfQ3Bins{"ConfQ3Bins", {2000, 0., 8.}, "binning Q3"};
  ConfigurableAxis ConfQ3BinsFor4D{"ConfQ3BinsFor4D", {500, 0., 2.}, "binning Q3 for 4D hist"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfFillCPRQA{"ConfFillCPRQA", false, "Fill Close Pair Rejection plots as a function of eta and phi"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};

  Configurable<float> ConfCPRdeltaPhiMax_pp{"ConfCPRdeltaPhiMax_pp", 0.01, "Max. Delta Phi for Close Pair Rejection of pp"};
  Configurable<float> ConfCPRdeltaEtaMax_pp{"ConfCPRdeltaEtaMax_pp", 0.01, "Max. Delta Eta for Close Pair Rejection of pp"};
  Configurable<float> ConfCPRdeltaPhiMax_pL{"ConfCPRdeltaPhiMax_pL", 0.1, "Max. Delta Phi for Close Pair Rejection of p and Lambda daughter"};
  Configurable<float> ConfCPRdeltaEtaMax_pL{"ConfCPRdeltaEtaMax_pL", 0.1, "Max. Delta Eta for Close Pair Rejection of p and Lambda daughter"};
  Configurable<float> ConfMaxQ3IncludedInCPRPlots{"ConfMaxQ3IncludedInCPRPlots", 8., "Maximum Q3, for which the pair CPR is included in plots"};
  ConfigurableAxis ConfDummy{"ConfDummy", {1, 0, 1}, "Dummy axis"};

  FemtoDreamContainerThreeBody<femtoDreamContainerThreeBody::EventType::same, femtoDreamContainerThreeBody::Observable::Q3> sameEventCont;
  FemtoDreamContainerThreeBody<femtoDreamContainerThreeBody::EventType::mixed, femtoDreamContainerThreeBody::Observable::Q3> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleanerTrackTrack;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCascade> pairCleanerTrackCascade;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejectionTrackTrackSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCascade> pairCloseRejectionTrackCascadeSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejectionTrackTrackME;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCascade> pairCloseRejectionTrackCascadeME;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry ThreeBodyQARegistry{"ThreeBodyQARegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext& context)
  {

    eventHisto.init(&qaRegistry, false);

    colBinning = {{ConfVtxBins, ConfMultBins}, true};

    trackHistoSelectedParts.init(&qaRegistry, ConfDummy, ConfDummy, ConfTempFitVarpTBins, ConfDummy, ConfDummy, ConfTempFitVarBinsTrack, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, ConfPDGCodePart);
    particleHistoSelectedCascades.init(&qaRegistry, ConfDummy, ConfDummy, ConfTempFitVarpTV0Bins, ConfDummy, ConfDummy, ConfTempFitVarBinsCascade, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfInvMassBins, ConfDummy, ConfIsMC, ConfPDGCodeCascade);
    particleHistoPosChild.init(&qaRegistry, ConfDummy, ConfDummy, ConfTempFitVarpTV0Child, ConfDummy, ConfDummy, ConfTempFitVarBinsCascadeChild, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, 0);
    particleHistoNegChild.init(&qaRegistry, ConfDummy, ConfDummy, ConfTempFitVarpTV0Child, ConfDummy, ConfDummy, ConfTempFitVarBinsCascadeChild, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, 0);
    particleHistoBachChild.init(&qaRegistry, ConfDummy, ConfDummy, ConfTempFitVarpTV0Child, ConfDummy, ConfDummy, ConfTempFitVarBinsCascadeChild, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, 0);

    ThreeBodyQARegistry.add("TripletTaskQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    ThreeBodyQARegistry.add("TripletTaskQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    ThreeBodyQARegistry.add("TripletTaskQA/hMinvSE_Cascade", ";Q_{3};M_{inv}", kTH2F, {ConfQ3Bins, ConfInvMassBins});
    ThreeBodyQARegistry.add("TripletTaskQA/hMinvME_Cascade", ";Q_{3};M_{inv}", kTH2F, {ConfQ3Bins, ConfInvMassBins});
    ThreeBodyQARegistry.add("TripletTaskQA/hMinvSE_LambdaChild", ";Q_{3};M_{inv}", kTH2F, {ConfQ3Bins, ConfInvMassBins});
    ThreeBodyQARegistry.add("TripletTaskQA/hMinvME_LambdaChild", ";Q_{3};M_{inv}", kTH2F, {ConfQ3Bins, ConfInvMassBins});
    ThreeBodyQARegistry.add("TripletTaskQA/particle_pT_in_Triplet_SE", "; p_{T1} ; p_{T2} ; p_{T3} ; Q_{3}", kTHnSparseF, {ConfTempFitVarpTBins, ConfTempFitVarpTBins, ConfTempFitVarpTV0Bins, ConfQ3BinsFor4D});
    ThreeBodyQARegistry.add("TripletTaskQA/particle_pT_in_Triplet_ME", "; p_{T1} ; p_{T2} ; p_{T3} ; Q_{3}", kTHnSparseF, {ConfTempFitVarpTBins, ConfTempFitVarpTBins, ConfTempFitVarpTV0Bins, ConfQ3BinsFor4D});
    ThreeBodyQARegistry.add("TripletTaskQA/hCentrality", ";Centrality; Q3", kTH2F, {{100, 0, 100}, ConfQ3Bins});

    std::vector<double> tmpVecMult = ConfMultBins;
    framework::AxisSpec multAxis = {tmpVecMult, "Multiplicity"};
    ThreeBodyQARegistry.add("TripletTaskQA/hSEMultVSGoodTracks", ";Mult;GoodT", kTH2F, {multAxis, {100, 0, 100}});
    ThreeBodyQARegistry.add("TripletTaskQA/hTestPairCleaner", ";posDaughtID; negDaughID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    ThreeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerPos", ";primaryTrack; posDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    ThreeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerNeg", ";primaryTrack; negDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    ThreeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerPosGlobal", ";primaryTrackGlobal; posDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    ThreeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerNegGlobal", ";primaryTrackGlobal; negDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    ThreeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerPosAfter", ";primaryTrack; posDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    ThreeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerNegAfter", ";primaryTrack; negDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});

    sameEventCont.init(&resultRegistry, ConfQ3Bins, ConfMultBins, ConfIsMC);
    mixedEventCont.init(&resultRegistry, ConfQ3Bins, ConfMultBins, ConfIsMC);
    sameEventCont.setPDGCodes(ConfPDGCodePart, ConfPDGCodePart, ConfPDGCodeCascade);
    mixedEventCont.setPDGCodes(ConfPDGCodePart, ConfPDGCodePart, ConfPDGCodeCascade);

    pairCleanerTrackTrack.init(&qaRegistry);
    pairCleanerTrackCascade.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejectionTrackTrackSE.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax_pp.value, ConfCPRdeltaEtaMax_pp.value, ConfCPRPlotPerRadii.value, 1, ConfUseOLD_possiblyWrong_CPR, ConfMaxQ3IncludedInCPRPlots, false, ConfAtWhichRadiiToCut, ConfAtWhichTPCRadii, ConfFillCPRQA);
      pairCloseRejectionTrackCascadeSE.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax_pL.value, ConfCPRdeltaEtaMax_pL.value, ConfCPRPlotPerRadii.value, 1, ConfUseOLD_possiblyWrong_CPR, ConfMaxQ3IncludedInCPRPlots, false, ConfAtWhichRadiiToCut, ConfAtWhichTPCRadii, ConfFillCPRQA);
      pairCloseRejectionTrackTrackME.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax_pp.value, ConfCPRdeltaEtaMax_pp.value, ConfCPRPlotPerRadii.value, 2, ConfUseOLD_possiblyWrong_CPR, ConfMaxQ3IncludedInCPRPlots, false, ConfAtWhichRadiiToCut, ConfAtWhichTPCRadii, ConfFillCPRQA);
      pairCloseRejectionTrackCascadeME.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax_pL.value, ConfCPRdeltaEtaMax_pL.value, ConfCPRPlotPerRadii.value, 2, ConfUseOLD_possiblyWrong_CPR, ConfMaxQ3IncludedInCPRPlots, true, ConfAtWhichRadiiToCut, ConfAtWhichTPCRadii, ConfFillCPRQA);
    }

    // get masses
    mMassOne = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePart)->Mass();
    mMassTwo = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePart)->Mass();
    mMassThree = TDatabasePDG::Instance()->GetParticle(ConfPDGCodeCascade)->Mass();

    // get bit for the collision mask
    std::bitset<8 * sizeof(aod::femtodreamcollision::BitMaskType)> mask;
    int index = 0;
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.find("femto-dream-triplet-task-track-track-cascade") != std::string::npos) {
        if (containsNameValuePair(device.options, "ConfCutPart", ConfCutPart.value) &&
            containsNameValuePair(device.options, "ConfTPCPIDBit", ConfTPCPIDBit.value) &&
            containsNameValuePair(device.options, "ConfTPCTOFPIDBit", ConfTPCTOFPIDBit.value) &&
            containsNameValuePair(device.options, "ConfPIDthrMom", ConfPIDthrMom.value) &&
            containsNameValuePair(device.options, "ConfMaxpT", ConfMaxpT.value) &&
            containsNameValuePair(device.options, "ConfMinpT", ConfMinpT.value) &&
            containsNameValuePair(device.options, "ConfCutCascade", ConfCutCascade.value) &&
            containsNameValuePair(device.options, "ChildPos_CutBit", ChildPos_CutBit.value) &&
            containsNameValuePair(device.options, "ChildPos_TPCBit", ChildPos_TPCBit.value) &&
            containsNameValuePair(device.options, "ChildNeg_CutBit", ChildNeg_CutBit.value) &&
            containsNameValuePair(device.options, "ChildNeg_TPCBit", ChildNeg_TPCBit.value) &&
            containsNameValuePair(device.options, "ChildBach_CutBit", ChildBach_CutBit.value) &&
            containsNameValuePair(device.options, "ChildBach_TPCBit", ChildBach_TPCBit.value) &&
            containsNameValuePair(device.options, "Conf_minInvMass_Cascade", Conf_minInvMass_Cascade.value) &&
            containsNameValuePair(device.options, "Conf_maxInvMass_Cascade", Conf_maxInvMass_Cascade.value) &&
            containsNameValuePair(device.options, "InvMassV0DaughMin", InvMassV0DaughMin.value) &&
            containsNameValuePair(device.options, "InvMassV0DaughMax", InvMassV0DaughMax.value) &&
            containsNameValuePair(device.options, "Conf_minPt_Cascade", Conf_minPt_Cascade.value) &&
            containsNameValuePair(device.options, "Conf_maxPt_Cascade", Conf_maxPt_Cascade.value)) {
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

    if ((ConfMixIfTripletPresent && ConfMixIfTVOPairPresent) ||
        (ConfMixIfTripletPresent && ConfMixIfTOrVOPartsPresent) ||
        (ConfMixIfTVOPairPresent && ConfMixIfTOrVOPartsPresent)) {
      LOG(fatal) << "Only one method of mixing can be chosen!";
    }
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
  void doSameEvent(PartitionType groupSelectedTracks, PartitionType groupSelectedCascades, PartType parts, float magFieldTesla, int multCol, float centCol)
  {
    /// Histograming tracks
    for (auto& part : groupSelectedTracks) {
      trackHistoSelectedParts.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt, multCol, centCol);
    }
    /// Histograming V0s
    for (auto& casc : groupSelectedCascades) {
      const auto& posChild = parts.iteratorAt(casc.index() - 3);
      const auto& negChild = parts.iteratorAt(casc.index() - 2);
      const auto& bachChild = parts.iteratorAt(casc.index() - 1);

      if (((posChild.cut() & ChildPos_CutBit) == ChildPos_CutBit &&
           (posChild.pidcut() & ChildPos_TPCBit) == ChildPos_TPCBit &&
           (negChild.cut() & ChildNeg_CutBit) == ChildNeg_CutBit &&
           (negChild.pidcut() & ChildNeg_TPCBit) == ChildNeg_TPCBit &&
           (bachChild.cut() & ChildBach_CutBit) == ChildBach_CutBit &&
           (bachChild.pidcut() & ChildBach_TPCBit) == ChildBach_TPCBit)) {
        particleHistoSelectedCascades.fillQA<isMC, false>(casc, aod::femtodreamparticle::kPt, multCol, centCol);
        particleHistoPosChild.fillQA<isMC, false>(posChild, aod::femtodreamparticle::kPt, multCol, centCol);
        particleHistoNegChild.fillQA<isMC, false>(negChild, aod::femtodreamparticle::kPt, multCol, centCol);
        particleHistoBachChild.fillQA<isMC, false>(bachChild, aod::femtodreamparticle::kPt, multCol, centCol);
      }
    }

    /// Now build the combinations
    for (auto& casc : groupSelectedCascades) {
      const auto& posChild = parts.iteratorAt(casc.index() - 3);
      const auto& negChild = parts.iteratorAt(casc.index() - 2);
      const auto& bachChild = parts.iteratorAt(casc.index() - 1);

      const auto& childrenPos = posChild.childrenIds();
      const auto& childrenNeg = negChild.childrenIds();
      const auto& childrenBac = bachChild.childrenIds();
      auto posID = childrenPos[0];
      auto negID = childrenNeg[1];
      auto bacID = childrenBac[2];
      ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleaner"), posID, negID);

      if (!((posChild.cut() & ChildPos_CutBit) == ChildPos_CutBit &&
            (posChild.pidcut() & ChildPos_TPCBit) == ChildPos_TPCBit &&
            (negChild.cut() & ChildNeg_CutBit) == ChildNeg_CutBit &&
            (negChild.pidcut() & ChildNeg_TPCBit) == ChildNeg_TPCBit)) {
        continue;
      }

      for (auto& [T1, T2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupSelectedTracks, groupSelectedTracks))) {
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPos"), T1.index(), posID);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNeg"), T1.index(), negID);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPos"), T2.index(), posID);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNeg"), T2.index(), negID);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPosGlobal"), T1.globalIndex(), posID);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNegGlobal"), T1.globalIndex(), negID);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPosGlobal"), T2.globalIndex(), posID);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNegGlobal"), T2.globalIndex(), negID);
        auto Q3 = FemtoDreamMath::getQ3(T1, mMassOne, T2, mMassTwo, casc, mMassThree);
        // Close pair rejection
        if (ConfIsCPR.value) {
          if (pairCloseRejectionTrackTrackSE.isClosePair(T1, T2, parts, magFieldTesla, Q3)) {
            continue;
          }
          if (pairCloseRejectionTrackCascadeSE.isClosePair(T1, casc, parts, magFieldTesla, Q3)) {
            continue;
          }
          if (pairCloseRejectionTrackCascadeSE.isClosePair(T2, casc, parts, magFieldTesla, Q3)) {
            continue;
          }
        }

        // track cleaning
        if (!pairCleanerTrackTrack.isCleanPair(T1, T2, parts)) {
          continue;
        }
        if (!pairCleanerTrackCascade.isCleanPair(T2, casc, parts)) {
          continue;
        }
        if (!pairCleanerTrackCascade.isCleanPair(T1, casc, parts)) {
          continue;
        }
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPosAfter"), T1.index(), posID);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNegAfter"), T1.index(), negID);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPosAfter"), T2.index(), posID);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNegAfter"), T2.index(), negID);
        // fill inv Mass as a function of Q3 for purity fits
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hMinvSE_Cascade"), Q3, casc.mLambda());
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hMinvSE_LambdaChild"), Q3, casc.mAntiLambda());
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/particle_pT_in_Triplet_SE"), T1.pt(), T2.pt(), casc.pt(), Q3);
        sameEventCont.setTriplet<isMC>(T1, T2, casc, multCol, Q3);
        ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hCentrality"), centCol, Q3);
      }
    }
  }

  /// process function to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEvent(o2::aod::FDCollision& col,
                        o2::aod::FDParticles& parts)
  {
    eventHisto.fillQA<false>(col);
    auto thegroupSelectedTracks = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedCascades = SelectedCascades->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache); /// Histograming V0s
    if (thegroupSelectedTracks.size() < 2 || thegroupSelectedCascades.size() < 1) {
      return;
    }
    doSameEvent<false>(thegroupSelectedTracks, thegroupSelectedCascades, parts, col.magField(), col.multNtr(), col.multV0M());
  }
  PROCESS_SWITCH(femtoDreamTripletTaskTrackTrackCascade, processSameEvent, "Enable processing same event", true);

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
    for (auto& [T1, T2, casc] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo, groupPartsThree))) {
      const auto& posChild = parts.iteratorAt(casc.globalIndex() - 3);
      const auto& negChild = parts.iteratorAt(casc.globalIndex() - 2);
      const auto& bachChild = parts.iteratorAt(casc.globalIndex() - 1);

      if (((posChild.cut() & ChildPos_CutBit) == ChildPos_CutBit &&
           (posChild.pidcut() & ChildPos_TPCBit) == ChildPos_TPCBit &&
           (negChild.cut() & ChildNeg_CutBit) == ChildNeg_CutBit &&
           (negChild.pidcut() & ChildNeg_TPCBit) == ChildNeg_TPCBit &&
           (bachChild.cut() & ChildBach_CutBit) == ChildBach_CutBit &&
           (bachChild.pidcut() & ChildBach_TPCBit) == ChildBach_TPCBit)) {
        continue;
      }

      auto Q3 = FemtoDreamMath::getQ3(T1, mMassOne, T2, mMassTwo, casc, mMassThree);
      // Close pair rejection
      if (ConfIsCPR.value) {
        if (pairCloseRejectionTrackTrackME.isClosePair(T1, T2, parts, magFieldTesla, Q3)) {
          continue;
        }
        if (pairCloseRejectionTrackCascadeME.isClosePair(T1, casc, parts, magFieldTesla, Q3)) {
          continue;
        }
        if (pairCloseRejectionTrackCascadeME.isClosePair(T2, casc, parts, magFieldTesla, Q3)) {
          continue;
        }
      }

      // fill inv Mass as a function of Q3 for purity fits
      ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hMinvME_Cascade"), Q3, casc.mLambda());
      ThreeBodyQARegistry.fill(HIST("TripletTaskQA/hMinvME_LambdaChild"), Q3, casc.mAntiLambda());
      ThreeBodyQARegistry.fill(HIST("TripletTaskQA/particle_pT_in_Triplet_ME"), T1.pt(), T2.pt(), casc.pt(), Q3);
      mixedEventCont.setTriplet<isMC>(T1, T2, casc, multCol, Q3);
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
      auto groupPartsThree = SelectedCascades->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();
      const auto& magFieldTesla3 = collision3.magField();

      if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
        continue;
      }

      doMixedEvent<false>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoDreamTripletTaskTrackTrackCascade, processMixedEvent, "Enable processing mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamTripletTaskTrackTrackCascade>(cfgc),
  };
  return workflow;
}
