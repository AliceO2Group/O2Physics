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

/// \file femtoDreamTripletTaskTrackTrackV0PbPb.cxx
/// \brief Tasks that reads the track and V0 tables and creates triplets; only two identical tracks and a V0 can be used
/// \author Laura Serksnyte, CERN, laura.serksnyte@cern.ch
/// \author Wioleta Rzęsa, TU München, wioleta.rzesa@cern.ch

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

#include <bitset>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct FemtoDreamTripletTaskTrackTrackV0PbPb {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  Configurable<float> confCentralityMin{"confCentralityMin", 0, "Event sel: Minimum Centrality Percentile"};
  Configurable<float> confCentralityMax{"confCentralityMax", 10, "Event sel: Maximum Centrality Percentile"};
  Configurable<float> confZVertexCut{"confZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};
  Filter eventCentrality = aod::femtodreamcollision::multV0M >= confCentralityMin && aod::femtodreamcollision::multV0M <= confCentralityMax;
  Filter eventVertex = (nabs(aod::collision::posZ) < confZVertexCut);
  using FilteredFDCollisions = soa::Filtered<aod::FDCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;
  using MaskedCollisions = soa::Filtered<soa::Join<aod::FDCollisions, aod::FDColMasks>>;
  using MaskedCollision = MaskedCollisions::iterator;

  aod::femtodreamcollision::BitMaskType maskBit = -1;
  float mMassOne = -999, mMassTwo = -999, mMassThree = -999;

  // Pair/triplet cuts
  Configurable<bool> confMixIfTripletPresent{"confMixIfTripletPresent", true, "Use for mixing only events which have a TTV0 triplet"};
  Configurable<bool> confMixIfTVOPairPresent{"confMixIfTVOPairPresent", false, "Use for mixing only events which have a TV0 pair (at least one track and one V0)"};
  Configurable<bool> confMixIfTOrVOPartsPresent{"confMixIfTOrVOPartsPresent", false, "Use for mixing only events which have at least one particle of interest"};
  Configurable<int> confMinTrackNumber{"confMinTrackNumber", 2, "Minimum number of tracks in the event"};
  Configurable<int> confMinV0Number{"confMinV0Number", 1, "Minimum number of V0 in the event"};

  // which CPR to use, old is with a possible bug and new is fixed
  Configurable<bool> confUseOLDPossiblyWrongCPR{"confUseOLDPossiblyWrongCPR", true, "Use for old CPR, which possibly has a bug. This is implemented only for debugging reasons to compare old and new code on hyperloop datasets."};
  Configurable<int> confAtWhichRadiiToCut{"confAtWhichRadiiToCut", 1, "At which radii perform deta dphi selection: 0 - at PV, 1 - averaged phi, 2 - at given radii"};
  Configurable<float> confAtWhichTPCRadii{"confAtWhichTPCRadii", 85., "If confAtWhichRadiiToCut = 2; this allows to select at which TPC radii to cut"};

  /// First 2 tracks of the triplet
  Configurable<int> confPDGCodePart{"confPDGCodePart", 2212, "Particle PDG code"}; // proton
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confCutPart{"confCutPart", 5542474, "Track - Selection bit from cutCulator"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confTPCPIDBit{"confTPCPIDBit", 16, "PID TPC bit from cutCulator "};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confTPCTOFPIDBit{"confTPCTOFPIDBit", 8, "PID TPCTOF bit from cutCulator"};
  Configurable<float> confMaxpT{"confMaxpT", 4.0f, "Maximum transverse momentum of the particles"};
  Configurable<float> confMinpT{"confMinpT", 0.5f, "Minimum transverse momentum of the particles"};
  Configurable<float> confMaxDCAxy{"confMaxDCAxy", 0.1f, "Maximum DCAxy of the particles"};
  Configurable<float> confMinDCAxy{"confMinDCAxy", -0.1f, "Minimum DCAxy of the particles"};
  Configurable<float> confPIDthrMom{"confPIDthrMom", 0.75f, "Momentum threshold from which TPC and TOF are required for PID"};
  Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<bool> confDCACutPtDep{"confDCACutPtDep", false, "Use pt dependent dca cut for tracks"};
  Configurable<float> confDCACutPtDepPar0{"confDCACutPtDepPar0", 0.0105, "Parameter par[0] of the pt dep cut, par[0] + par[1]/(pT/(GeV/c)−1.1) cm"};
  Configurable<float> confDCACutPtDepPar1{"confDCACutPtDepPar1", 0.035, "Parameter par[1] of the pt dep cut, par[0] + par[1]/(pT/(GeV/c)−1.1) cm"};

  /// Partition for selected particles
  Partition<aod::FDParticles> selectedParts = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                              ifnode(aod::femtodreamparticle::pt * (0.5f * nexp(aod::femtodreamparticle::eta) + 0.5f * nexp(-1.f * aod::femtodreamparticle::eta)) <= confPIDthrMom, ncheckbit(aod::femtodreamparticle::pidcut, confTPCPIDBit), ncheckbit(aod::femtodreamparticle::pidcut, confTPCTOFPIDBit)) &&
                                              (ncheckbit(aod::femtodreamparticle::cut, confCutPart)) &&
                                              (aod::femtodreamparticle::pt < confMaxpT) &&
                                              (aod::femtodreamparticle::pt > confMinpT) &&
                                              ifnode(confDCACutPtDep, (nabs(aod::femtodreamparticle::tempFitVar) <= confDCACutPtDepPar0 + (confDCACutPtDepPar1 / npow(aod::femtodreamparticle::pt, 1.1f))),
                                                     ((aod::femtodreamparticle::tempFitVar >= confMinDCAxy) &&
                                                      (aod::femtodreamparticle::tempFitVar <= confMaxDCAxy)));
  ;

  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> selectedPartsMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                                            ifnode(aod::femtodreamparticle::pt * (0.5f * nexp(aod::femtodreamparticle::eta) + 0.5f * nexp(-1.f * aod::femtodreamparticle::eta)) <= confPIDthrMom, ncheckbit(aod::femtodreamparticle::pidcut, confTPCPIDBit), ncheckbit(aod::femtodreamparticle::pidcut, confTPCTOFPIDBit)) &&
                                                                            (ncheckbit(aod::femtodreamparticle::cut, confCutPart)) &&
                                                                            (aod::femtodreamparticle::pt < confMaxpT) &&
                                                                            (aod::femtodreamparticle::pt > confMinpT) &&
                                                                            ifnode(confDCACutPtDep, (nabs(aod::femtodreamparticle::tempFitVar) <= confDCACutPtDepPar0 + (confDCACutPtDepPar1 / npow(aod::femtodreamparticle::pt, 1.1f))),
                                                                                   ((aod::femtodreamparticle::tempFitVar >= confMinDCAxy) &&
                                                                                    (aod::femtodreamparticle::tempFitVar <= confMaxDCAxy)));
  ;

  /// Histogramming of selected tracks
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoselectedParts;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 5> trackHistoALLselectedParts;

  /// V0 selection
  Configurable<int> confPDGCodeV0{"confPDGCodeV0", 3122, "V0 PDG code"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confCutV0{"confCutV0", 7518, "V0 - Selection bit from cutCulator"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confChildPosCutV0{"confChildPosCutV0", 8234, "Selection bit for positive child of V0"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confChildPosTPCBitV0{"confChildPosTPCBitV0", 1024, "PID TPC bit for positive child of V0"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confChildNegCutV0{"confChildNegCutV0", 8233, "Selection bit for negative child of V0"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confChildNegTPCBitV0{"confChildNegTPCBitV0", 4096, "PID TPC bit for negative child of V0"};

  Configurable<float> confMinInvMassV0{"confMinInvMassV0", 1.08, "Minimum invariant mass of V0 (particle)"};
  Configurable<float> confMaxInvMassV0{"confMaxInvMassV0", 1.15, "Maximum invariant mass of V0 (particle)"};
  Configurable<float> confMinInvMassAntiV0{"confMinInvMassAntiV0", 1.08, "Minimum invariant mass of V0 (antiparticle)"};
  Configurable<float> confMaxInvMassAntiV0{"confMaxInvMassAntiV0", 1.15, "Maximum invariant mass of V0 (antiparticle)"};

  Configurable<float> confMinPtV0{"confMinPtV0", 0.0, "Minimum pT of V0"};
  Configurable<float> confMaxPtV0{"confMaxPtV0", 999.0, "Maximum pT of V0"};

  // Partition for selected particles
  Partition<aod::FDParticles> selectedV0s = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                            (ncheckbit(aod::femtodreamparticle::cut, confCutV0)) &&
                                            (aod::femtodreamparticle::mLambda > confMinInvMassV0) &&
                                            (aod::femtodreamparticle::mLambda < confMaxInvMassV0) &&
                                            (aod::femtodreamparticle::mAntiLambda > confMinInvMassAntiV0) &&
                                            (aod::femtodreamparticle::mAntiLambda < confMaxInvMassAntiV0) &&
                                            (aod::femtodreamparticle::pt > confMinPtV0) &&
                                            (aod::femtodreamparticle::pt < confMaxPtV0);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> selectedV0sMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                                                          (ncheckbit(aod::femtodreamparticle::cut, confCutV0)) &&
                                                                          (aod::femtodreamparticle::mLambda > confMinInvMassV0) &&
                                                                          (aod::femtodreamparticle::mLambda < confMaxInvMassV0) &&
                                                                          (aod::femtodreamparticle::mAntiLambda > confMinInvMassAntiV0) &&
                                                                          (aod::femtodreamparticle::mAntiLambda < confMaxInvMassAntiV0) &&
                                                                          (aod::femtodreamparticle::pt > confMinPtV0) &&
                                                                          (aod::femtodreamparticle::pt < confMaxPtV0);

  /// Histogramming of selected V0s
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 2> particleHistoselectedV0s;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 3> particleHistoPosChild;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 4> particleHistoNegChild;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 5> particleHistoALLselectedV0s;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 6> particleHistoALLPosChild;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 7> particleHistoALLNegChild;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// Particles
  ConfigurableAxis confTempFitVarBinsTrack{"confTempFitVarBinsTrack", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarpTBins{"confTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (track)"};
  ConfigurableAxis confTempFitVarBinsV0{"confTempFitVarBinsV0", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis confTempFitVarBinsV0Child{"confTempFitVarBinsV0Child", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0 child)"};
  ConfigurableAxis confTempFitVarpTV0Bins{"confTempFitVarpTV0Bins", {35, 0, 6}, "pT binning of the pT vs. TempFitVar plot (V0)"};
  ConfigurableAxis confTempFitVarpTV0Child{"confTempFitVarpTV0Child", {35, 0, 6}, "pT binning of the pT vs. TempFitVar plot (V0 child)"};
  ConfigurableAxis confInvMassBins{"confInvMassBins", {200, 1, 1.2}, "InvMass binning"};

  /// Correlations
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 23.0f, 38.0f, 53.0f, 81.0f, 110.0f, 157.0f, 205.0f, 278.0f, 351.0f, 455.0f, 559.0f, 703.0f, 848.0f, 1050.0f, 1253.0f, 1530.0f, 1668.0f, 1857.0f, 2047.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

  ConfigurableAxis confQ3Bins{"confQ3Bins", {2000, 0., 8.}, "binning Q3"};
  ConfigurableAxis confQ3BinsFor4D{"confQ3BinsFor4D", {500, 0., 2.}, "binning Q3 for 4D hist"};
  Configurable<int> confNEventsMix{"confNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> confIsCPR{"confIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> confFillCPRQA{"confFillCPRQA", false, "Fill Close Pair Rejection plots as a function of eta and phi"};
  Configurable<bool> confCPRPlotPerRadii{"confCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> confCPRdeltaPhiMaxpp{"confCPRdeltaPhiMaxpp", 0.01, "Max. Delta Phi for Close Pair Rejection of pp"};
  Configurable<float> confCPRdeltaEtaMaxpp{"confCPRdeltaEtaMaxpp", 0.01, "Max. Delta Eta for Close Pair Rejection of pp"};
  Configurable<float> confCPRdeltaPhiMaxpL{"confCPRdeltaPhiMaxpL", 0.1, "Max. Delta Phi for Close Pair Rejection of p and Lambda daughter"};
  Configurable<float> confCPRdeltaEtaMaxpL{"confCPRdeltaEtaMaxpL", 0.1, "Max. Delta Eta for Close Pair Rejection of p and Lambda daughter"};
  Configurable<float> confMaxQ3IncludedInCPRPlots{"confMaxQ3IncludedInCPRPlots", 8., "Maximum Q3, for which the pair CPR is included in plots"};
  ConfigurableAxis confDummy{"confDummy", {1, 0, 1}, "Dummy axis"};

  FemtoDreamContainerThreeBody<femtoDreamContainerThreeBody::EventType::same, femtoDreamContainerThreeBody::Observable::Q3> sameEventCont;
  FemtoDreamContainerThreeBody<femtoDreamContainerThreeBody::EventType::mixed, femtoDreamContainerThreeBody::Observable::Q3> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleanerTrackTrack;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCleanerTrackV0;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejectionTrackTrackSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCloseRejectionTrackV0SE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejectionTrackTrackME;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCloseRejectionTrackV0ME;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry threeBodyQARegistry{"threeBodyQARegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext& context)
  {

    eventHisto.init(&qaRegistry, false);

    colBinning = {{confVtxBins, confMultBins}, true};

    trackHistoselectedParts.init(&qaRegistry, confDummy, confDummy, confTempFitVarpTBins, confDummy, confDummy, confTempFitVarBinsTrack, confDummy, confDummy, confDummy, confDummy, confDummy, confDummy, confIsMC, confPDGCodePart);
    trackHistoALLselectedParts.init(&qaRegistry, confDummy, confDummy, confTempFitVarpTBins, confDummy, confDummy, confTempFitVarBinsTrack, confDummy, confDummy, confDummy, confDummy, confDummy, confDummy, confIsMC, confPDGCodePart);
    particleHistoselectedV0s.init(&qaRegistry, confDummy, confDummy, confTempFitVarpTV0Bins, confDummy, confDummy, confTempFitVarBinsV0, confDummy, confDummy, confDummy, confDummy, confInvMassBins, confDummy, confIsMC, confPDGCodeV0);
    particleHistoPosChild.init(&qaRegistry, confDummy, confDummy, confTempFitVarpTV0Child, confDummy, confDummy, confTempFitVarBinsV0Child, confDummy, confDummy, confDummy, confDummy, confDummy, confDummy, confIsMC, 0);
    particleHistoNegChild.init(&qaRegistry, confDummy, confDummy, confTempFitVarpTV0Child, confDummy, confDummy, confTempFitVarBinsV0Child, confDummy, confDummy, confDummy, confDummy, confDummy, confDummy, confIsMC, 0);
    particleHistoALLselectedV0s.init(&qaRegistry, confDummy, confDummy, confTempFitVarpTV0Bins, confDummy, confDummy, confTempFitVarBinsV0, confDummy, confDummy, confDummy, confDummy, confInvMassBins, confDummy, confIsMC, confPDGCodeV0);
    particleHistoALLPosChild.init(&qaRegistry, confDummy, confDummy, confTempFitVarpTV0Child, confDummy, confDummy, confTempFitVarBinsV0Child, confDummy, confDummy, confDummy, confDummy, confDummy, confDummy, confIsMC, 0);
    particleHistoALLNegChild.init(&qaRegistry, confDummy, confDummy, confTempFitVarpTV0Child, confDummy, confDummy, confTempFitVarBinsV0Child, confDummy, confDummy, confDummy, confDummy, confDummy, confDummy, confIsMC, 0);

    threeBodyQARegistry.add("TripletTaskQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    threeBodyQARegistry.add("TripletTaskQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    threeBodyQARegistry.add("TripletTaskQA/hMinvSE_Lambda", ";Q_{3};M_{inv}", kTH2F, {confQ3Bins, confInvMassBins});
    threeBodyQARegistry.add("TripletTaskQA/hMinvME_Lambda", ";Q_{3};M_{inv}", kTH2F, {confQ3Bins, confInvMassBins});
    threeBodyQARegistry.add("TripletTaskQA/hMinvSE_AntiLambda", ";Q_{3};M_{inv}", kTH2F, {confQ3Bins, confInvMassBins});
    threeBodyQARegistry.add("TripletTaskQA/hMinvME_AntiLambda", ";Q_{3};M_{inv}", kTH2F, {confQ3Bins, confInvMassBins});
    threeBodyQARegistry.add("TripletTaskQA/particle_pT_in_Triplet_SE", "; p_{T1} ; p_{T2} ; p_{T3} ; Q_{3}", kTHnSparseF, {confTempFitVarpTBins, confTempFitVarpTBins, confTempFitVarpTV0Bins, confQ3BinsFor4D});
    threeBodyQARegistry.add("TripletTaskQA/particle_pT_in_Triplet_ME", "; p_{T1} ; p_{T2} ; p_{T3} ; Q_{3}", kTHnSparseF, {confTempFitVarpTBins, confTempFitVarpTBins, confTempFitVarpTV0Bins, confQ3BinsFor4D});
    threeBodyQARegistry.add("TripletTaskQA/hCentrality", ";Centrality; Q3", kTH2F, {{100, 0, 100}, confQ3Bins});

    std::vector<double> tmpVecMult = confMultBins;
    framework::AxisSpec multAxis = {tmpVecMult, "Multiplicity"};
    threeBodyQARegistry.add("TripletTaskQA/hSEMultVSGoodTracks", ";Mult;GoodT", kTH2F, {multAxis, {100, 0, 100}});
    threeBodyQARegistry.add("TripletTaskQA/hTestPairCleaner", ";posDaughtID; negDaughID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    threeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerPos", ";primaryTrack; posDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    threeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerNeg", ";primaryTrack; negDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    threeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerPosGlobal", ";primaryTrackGlobal; posDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    threeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerNegGlobal", ";primaryTrackGlobal; negDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    threeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerPosAfter", ";primaryTrack; posDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    threeBodyQARegistry.add("TripletTaskQA/hTestPairCleanerNegAfter", ";primaryTrack; negDaughtID", kTH2F, {{40, -20, 20}, {40, -20, 20}});
    threeBodyQARegistry.add("TripletTaskQA/hCentralityME", ";Centrality;Entries", kTH1F, {{100, 0.0, 100.0}});

    threeBodyQARegistry.add("SameEvent/relPairDist_trackTrack", ";k* (GeV/c) ;Entries", kTH1F, {{2000, 0.0, 4.0}});
    threeBodyQARegistry.add("MixedEvent/relPairDist_trackTrack", ";k* (GeV/c) ;Entries", kTH1F, {{2000, 0.0, 4.0}});
    threeBodyQARegistry.add("SameEvent/relPairDist_track1V0", ";k* (GeV/c) ;Entries", kTH1F, {{2000, 0.0, 4.0}});
    threeBodyQARegistry.add("SameEvent/relPairDist_track2V0", ";k* (GeV/c) ;Entries", kTH1F, {{2000, 0.0, 4.0}});
    threeBodyQARegistry.add("MixedEvent/relPairDist_track1V0", ";k* (GeV/c) ;Entries", kTH1F, {{2000, 0.0, 4.0}});
    threeBodyQARegistry.add("MixedEvent/relPairDist_track2V0", ";k* (GeV/c) ;Entries", kTH1F, {{2000, 0.0, 4.0}});
    sameEventCont.init(&resultRegistry, confQ3Bins, confMultBins, confIsMC);
    mixedEventCont.init(&resultRegistry, confQ3Bins, confMultBins, confIsMC);
    sameEventCont.setPDGCodes(confPDGCodePart, confPDGCodePart, confPDGCodeV0);
    mixedEventCont.setPDGCodes(confPDGCodePart, confPDGCodePart, confPDGCodeV0);

    pairCleanerTrackTrack.init(&qaRegistry);
    pairCleanerTrackV0.init(&qaRegistry);
    if (confIsCPR.value) {
      pairCloseRejectionTrackTrackSE.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiMaxpp.value, confCPRdeltaEtaMaxpp.value, confCPRPlotPerRadii.value, 1, confUseOLDPossiblyWrongCPR, confMaxQ3IncludedInCPRPlots, false, confAtWhichRadiiToCut, confAtWhichTPCRadii, confFillCPRQA);
      pairCloseRejectionTrackV0SE.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiMaxpL.value, confCPRdeltaEtaMaxpL.value, confCPRPlotPerRadii.value, 1, confUseOLDPossiblyWrongCPR, confMaxQ3IncludedInCPRPlots, false, confAtWhichRadiiToCut, confAtWhichTPCRadii, confFillCPRQA);
      pairCloseRejectionTrackTrackME.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiMaxpp.value, confCPRdeltaEtaMaxpp.value, confCPRPlotPerRadii.value, 2, confUseOLDPossiblyWrongCPR, confMaxQ3IncludedInCPRPlots, false, confAtWhichRadiiToCut, confAtWhichTPCRadii, confFillCPRQA);
      pairCloseRejectionTrackV0ME.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiMaxpL.value, confCPRdeltaEtaMaxpL.value, confCPRPlotPerRadii.value, 2, confUseOLDPossiblyWrongCPR, confMaxQ3IncludedInCPRPlots, true, confAtWhichRadiiToCut, confAtWhichTPCRadii, confFillCPRQA);
    }

    // get masses
    mMassOne = o2::constants::physics::MassProton;
    mMassTwo = o2::constants::physics::MassProton;
    mMassThree = o2::constants::physics::MassLambda;

    // get bit for the collision mask
    std::bitset<8 * sizeof(aod::femtodreamcollision::BitMaskType)> mask;
    int index = 0;
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.find("femto-dream-triplet-task-track-track-v0-pb-pb") != std::string::npos) {
        if (containsNameValuePair(device.options, "confCutPart", confCutPart.value) &&
            containsNameValuePair(device.options, "confTPCPIDBit", confTPCPIDBit.value) &&
            containsNameValuePair(device.options, "confTPCTOFPIDBit", confTPCTOFPIDBit.value) &&
            containsNameValuePair(device.options, "confPIDthrMom", confPIDthrMom.value) &&
            containsNameValuePair(device.options, "confMaxpT", confMaxpT.value) &&
            containsNameValuePair(device.options, "confMinpT", confMinpT.value) &&
            containsNameValuePair(device.options, "confCutV0", confCutV0.value) &&
            containsNameValuePair(device.options, "confChildPosCutV0", confChildPosCutV0.value) &&
            containsNameValuePair(device.options, "confChildPosTPCBitV0", confChildPosTPCBitV0.value) &&
            containsNameValuePair(device.options, "confChildNegCutV0", confChildNegCutV0.value) &&
            containsNameValuePair(device.options, "confChildNegTPCBitV0", confChildNegTPCBitV0.value) &&
            containsNameValuePair(device.options, "confMinInvMassV0", confMinInvMassV0.value) &&
            containsNameValuePair(device.options, "confMaxInvMassV0", confMaxInvMassV0.value) &&
            containsNameValuePair(device.options, "confMinInvMassAntiV0", confMinInvMassAntiV0.value) &&
            containsNameValuePair(device.options, "confMaxInvMassAntiV0", confMaxInvMassAntiV0.value) &&
            containsNameValuePair(device.options, "confMinPtV0", confMinPtV0.value) &&
            containsNameValuePair(device.options, "confMaxPtV0", confMaxPtV0.value)) {
          mask.set(index);
          maskBit = static_cast<aod::femtodreamcollision::BitMaskType>(mask.to_ulong());
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

    if ((confMixIfTripletPresent && confMixIfTVOPairPresent) ||
        (confMixIfTripletPresent && confMixIfTOrVOPartsPresent) ||
        (confMixIfTVOPairPresent && confMixIfTOrVOPartsPresent)) {
      LOG(fatal) << "Only one method of mixing can be chosen!";
    }
  }

  template <bool isMC, typename CollisionType>
  void fillCollision(CollisionType col)
  {
    threeBodyQARegistry.fill(HIST("TripletTaskQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
    eventHisto.fillQA<isMC>(col);
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
  void doSameEvent(PartitionType groupSelectedTracks, PartitionType groupselectedV0s, PartType parts, float magFieldTesla, int multCol, float centCol)
  {
    /// Histograming tracks
    for (const auto& part : groupSelectedTracks) {
      trackHistoselectedParts.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt, multCol, centCol);
    }
    /// Histograming V0s
    for (const auto& V0 : groupselectedV0s) {
      const auto& posChild = parts.iteratorAt(V0.index() - 2);
      const auto& negChild = parts.iteratorAt(V0.index() - 1);

      if (((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
           (posChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0 &&
           (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
           (negChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0)) {
        particleHistoselectedV0s.fillQA<isMC, false>(V0, aod::femtodreamparticle::kPt, multCol, centCol);
        particleHistoPosChild.fillQA<isMC, false>(posChild, aod::femtodreamparticle::kPt, multCol, centCol);
        particleHistoNegChild.fillQA<isMC, false>(negChild, aod::femtodreamparticle::kPt, multCol, centCol);
      }
    }

    /// Now build the combinations
    for (const auto& V0 : groupselectedV0s) {
      const auto& posChild = parts.iteratorAt(V0.index() - 2);
      const auto& negChild = parts.iteratorAt(V0.index() - 1);

      const auto& childrenPos = posChild.childrenIds();
      const auto& childrenNeg = negChild.childrenIds();
      auto posID = childrenPos[0];
      auto negID = childrenNeg[1];
      threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleaner"), posID, negID);

      if (!((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
            (posChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0 &&
            (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
            (negChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0)) {
        continue;
      }

      for (const auto& [T1, T2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupSelectedTracks, groupSelectedTracks))) {
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPos"), T1.index(), posID);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNeg"), T1.index(), negID);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPos"), T2.index(), posID);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNeg"), T2.index(), negID);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPosGlobal"), T1.globalIndex(), posID);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNegGlobal"), T1.globalIndex(), negID);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPosGlobal"), T2.globalIndex(), posID);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNegGlobal"), T2.globalIndex(), negID);
        auto q3 = FemtoDreamMath::getQ3(T1, mMassOne, T2, mMassTwo, V0, mMassThree);
        // Close pair rejection
        if (confIsCPR.value) {
          if (pairCloseRejectionTrackTrackSE.isClosePair(T1, T2, parts, magFieldTesla, q3)) {
            continue;
          }
          if (pairCloseRejectionTrackV0SE.isClosePair(T1, V0, parts, magFieldTesla, q3)) {
            continue;
          }
          if (pairCloseRejectionTrackV0SE.isClosePair(T2, V0, parts, magFieldTesla, q3)) {
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
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPosAfter"), T1.index(), posID);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNegAfter"), T1.index(), negID);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerPosAfter"), T2.index(), posID);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hTestPairCleanerNegAfter"), T2.index(), negID);
        // fill inv Mass as a function of Q3 for purity fits
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hMinvSE_Lambda"), q3, V0.mLambda());
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hMinvSE_AntiLambda"), q3, V0.mAntiLambda());
        threeBodyQARegistry.fill(HIST("TripletTaskQA/particle_pT_in_Triplet_SE"), T1.pt(), T2.pt(), V0.pt(), q3);
        sameEventCont.setTriplet<isMC>(T1, T2, V0, multCol, q3);
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hCentrality"), centCol, q3);

        auto kstarTT = FemtoDreamMath::getkstar(T1, mMassOne, T2, mMassTwo);
        auto kstarT1V = FemtoDreamMath::getkstar(T1, mMassOne, V0, mMassThree);
        auto kstarT2V = FemtoDreamMath::getkstar(T2, mMassTwo, V0, mMassThree);
        threeBodyQARegistry.fill(HIST("SameEvent/relPairDist_trackTrack"), kstarTT);
        threeBodyQARegistry.fill(HIST("SameEvent/relPairDist_track1V0"), kstarT1V);
        threeBodyQARegistry.fill(HIST("SameEvent/relPairDist_track2V0"), kstarT2V);
      }
    }
  }

  /// process function to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEvent(const FilteredFDCollision& col, const o2::aod::FDParticles& parts)
  {
    fillCollision<false>(col);
    auto thegroupSelectedTracks = selectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    for (const auto& part : thegroupSelectedTracks) {
      trackHistoALLselectedParts.fillQA<false, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
    }
    auto thegroupselectedV0s = selectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache); /// Histograming V0s
    for (const auto& V0 : thegroupselectedV0s) {
      const auto& posChild = parts.iteratorAt(V0.index() - 2);
      const auto& negChild = parts.iteratorAt(V0.index() - 1);

      if (((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
           (posChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0 &&
           (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
           (negChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0)) {
        particleHistoALLselectedV0s.fillQA<false, false>(V0, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        particleHistoALLPosChild.fillQA<false, false>(posChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        particleHistoALLNegChild.fillQA<false, false>(negChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }
    if (thegroupSelectedTracks.size() < confMinTrackNumber || thegroupselectedV0s.size() < confMinV0Number) {
      return;
    }
    doSameEvent<false>(thegroupSelectedTracks, thegroupselectedV0s, parts, col.magField(), col.multNtr(), col.multV0M());
  }
  PROCESS_SWITCH(FemtoDreamTripletTaskTrackTrackV0PbPb, processSameEvent, "Enable processing same event", true);

  /// process function to call doSameEvent with Data which has a mask for containing particles or not
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEventMasked(const MaskedCollision& col, const o2::aod::FDParticles& parts)
  {
    fillCollision<false>(col);
    auto thegroupSelectedTracks = selectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    for (const auto& part : thegroupSelectedTracks) {
      trackHistoALLselectedParts.fillQA<false, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
    }
    auto thegroupselectedV0s = selectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    for (const auto& V0 : thegroupselectedV0s) {
      const auto& posChild = parts.iteratorAt(V0.index() - 2);
      const auto& negChild = parts.iteratorAt(V0.index() - 1);

      if (((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
           (posChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0 &&
           (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
           (negChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0)) {
        particleHistoALLselectedV0s.fillQA<false, false>(V0, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        particleHistoALLPosChild.fillQA<false, false>(posChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        particleHistoALLNegChild.fillQA<false, false>(negChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }
    if (thegroupSelectedTracks.size() < confMinTrackNumber || thegroupselectedV0s.size() < confMinV0Number) {
      return;
    }
    doSameEvent<false>(thegroupSelectedTracks, thegroupselectedV0s, parts, col.magField(), col.multNtr(), col.multV0M());
  }
  PROCESS_SWITCH(FemtoDreamTripletTaskTrackTrackV0PbPb, processSameEventMasked, "Enable processing same event with masks", false);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(const o2::aod::FDCollision& col,
                          const soa::Join<o2::aod::FDParticles, const o2::aod::FDMCLabels>& parts,
                          const o2::aod::FDMCParticles&)
  {
    fillCollision<false>(col);
    auto thegroupSelectedTracks = selectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    for (const auto& part : thegroupSelectedTracks) {
      trackHistoALLselectedParts.fillQA<true, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      threeBodyQARegistry.fill(HIST("TrackMC_QA/hMazzachi"), part.fdMCParticle().pt(), (part.pt() - part.fdMCParticle().pt()) / part.fdMCParticle().pt());
    }
    auto thegroupselectedV0s = selectedV0sMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    for (const auto& V0 : thegroupselectedV0s) {
      const auto& posChild = parts.iteratorAt(V0.index() - 2);
      const auto& negChild = parts.iteratorAt(V0.index() - 1);

      if (((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
           (posChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0 &&
           (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
           (negChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0)) {
        particleHistoALLselectedV0s.fillQA<true, false>(V0, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        particleHistoALLPosChild.fillQA<true, false>(posChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        particleHistoALLNegChild.fillQA<true, false>(negChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }
    if (thegroupSelectedTracks.size() < confMinTrackNumber || thegroupselectedV0s.size() < confMinV0Number) {
      return;
    }
    doSameEvent<true>(thegroupSelectedTracks, thegroupselectedV0s, parts, col.magField(), col.multNtr(), col.multV0M());
  }
  PROCESS_SWITCH(FemtoDreamTripletTaskTrackTrackV0PbPb, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  /// process function for to call doSameEvent with Monte Carlo which has a mask for containing particles or not
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMCMasked(const MaskedCollision& col,
                                const soa::Join<o2::aod::FDParticles, const o2::aod::FDMCLabels>& parts,
                                const o2::aod::FDMCParticles&)
  {
    fillCollision<false>(col);
    auto thegroupSelectedTracks = selectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    for (const auto& part : thegroupSelectedTracks) {
      trackHistoALLselectedParts.fillQA<true, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      threeBodyQARegistry.fill(HIST("TrackMC_QA/hMazzachi"), part.fdMCParticle().pt(), (part.pt() - part.fdMCParticle().pt()) / part.fdMCParticle().pt());
    }
    auto thegroupselectedV0s = selectedV0sMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    for (const auto& V0 : thegroupselectedV0s) {
      const auto& posChild = parts.iteratorAt(V0.index() - 2);
      const auto& negChild = parts.iteratorAt(V0.index() - 1);

      if (((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
           (posChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0 &&
           (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
           (negChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0)) {
        particleHistoALLselectedV0s.fillQA<true, false>(V0, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        particleHistoALLPosChild.fillQA<true, false>(posChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        particleHistoALLNegChild.fillQA<true, false>(negChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }
    if (thegroupSelectedTracks.size() < confMinTrackNumber || thegroupselectedV0s.size() < confMinV0Number) {
      return;
    }
    doSameEvent<true>(thegroupSelectedTracks, thegroupselectedV0s, parts, col.magField(), col.multNtr(), col.multV0M());
  }
  PROCESS_SWITCH(FemtoDreamTripletTaskTrackTrackV0PbPb, processSameEventMCMasked, "Enable processing same event for Monte Carlo", false);

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
    for (const auto& [T1, T2, V0] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo, groupPartsThree))) {
      const auto& posChild = parts.iteratorAt(V0.globalIndex() - 2);
      const auto& negChild = parts.iteratorAt(V0.globalIndex() - 1);

      if (!((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
            (posChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0 &&
            (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
            (negChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0)) {
        continue;
      }

      auto q3 = FemtoDreamMath::getQ3(T1, mMassOne, T2, mMassTwo, V0, mMassThree);
      // Close pair rejection
      if (confIsCPR.value) {
        if (pairCloseRejectionTrackTrackME.isClosePair(T1, T2, parts, magFieldTesla, q3)) {
          continue;
        }
        if (pairCloseRejectionTrackV0ME.isClosePair(T1, V0, parts, magFieldTesla, q3)) {
          continue;
        }
        if (pairCloseRejectionTrackV0ME.isClosePair(T2, V0, parts, magFieldTesla, q3)) {
          continue;
        }
      }

      // fill inv Mass as a function of Q3 for purity fits
      threeBodyQARegistry.fill(HIST("TripletTaskQA/hMinvME_Lambda"), q3, V0.mLambda());
      threeBodyQARegistry.fill(HIST("TripletTaskQA/hMinvME_AntiLambda"), q3, V0.mAntiLambda());
      threeBodyQARegistry.fill(HIST("TripletTaskQA/particle_pT_in_Triplet_ME"), T1.pt(), T2.pt(), V0.pt(), q3);
      mixedEventCont.setTriplet<isMC>(T1, T2, V0, multCol, q3);

      auto kstarTT = FemtoDreamMath::getkstar(T1, mMassOne, T2, mMassTwo);
      auto kstarT1V = FemtoDreamMath::getkstar(T1, mMassOne, V0, mMassThree);
      auto kstarT2V = FemtoDreamMath::getkstar(T2, mMassTwo, V0, mMassThree);
      threeBodyQARegistry.fill(HIST("MixedEvent/relPairDist_trackTrack"), kstarTT);
      threeBodyQARegistry.fill(HIST("MixedEvent/relPairDist_track1V0"), kstarT1V);
      threeBodyQARegistry.fill(HIST("MixedEvent/relPairDist_track2V0"), kstarT2V);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoDreamParticleTable
  void processMixedEvent(const FilteredFDCollisions& cols,
                         const o2::aod::FDParticles& parts)
  {
    for (const auto& [collision1, collision2, collision3] : soa::selfCombinations(colBinning, confNEventsMix, -1, cols, cols, cols)) {
      const int multiplicityCol = collision1.multNtr();
      threeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));
      threeBodyQARegistry.fill(HIST("TripletTaskQA/hCentralityME"), collision1.multV0M());
      threeBodyQARegistry.fill(HIST("TripletTaskQA/hCentralityME"), collision2.multV0M());
      threeBodyQARegistry.fill(HIST("TripletTaskQA/hCentralityME"), collision3.multV0M());

      auto groupPartsOne = selectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = selectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsThree = selectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();
      const auto& magFieldTesla3 = collision3.magField();
      if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
        continue;
      }

      doMixedEvent<false>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(FemtoDreamTripletTaskTrackTrackV0PbPb, processMixedEvent, "Enable processing mixed events", true);

  /// process function for to call doMixedEvent with Data which has a mask for containing particles or not
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoDreamParticleTable
  void processMixedEventMasked(const MaskedCollisions& cols, const o2::aod::FDParticles& parts)
  {
    if (confMixIfTripletPresent || confMixIfTVOPairPresent) {
      Partition<MaskedCollisions> partitionMaskedCol1 = (confMixIfTripletPresent && (aod::femtodreamcollision::bitmaskTrackTwo & maskBit) == maskBit && (aod::femtodreamcollision::bitmaskTrackThree & maskBit) == maskBit) ||
                                                        (confMixIfTVOPairPresent && (aod::femtodreamcollision::bitmaskTrackOne & maskBit) == maskBit && (aod::femtodreamcollision::bitmaskTrackThree & maskBit) == maskBit);
      partitionMaskedCol1.bindTable(cols);

      for (const auto& [collision1, collision2, collision3] : soa::selfCombinations(colBinning, confNEventsMix, -1, *partitionMaskedCol1.mFiltered, *partitionMaskedCol1.mFiltered, *partitionMaskedCol1.mFiltered)) {

        const int multiplicityCol = collision1.multNtr();
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hCentralityME"), collision1.multV0M());
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hCentralityME"), collision2.multV0M());
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hCentralityME"), collision3.multV0M());

        auto groupPartsOne = selectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = selectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
        auto groupPartsThree = selectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

        const auto& magFieldTesla1 = collision1.magField();
        const auto& magFieldTesla2 = collision2.magField();
        const auto& magFieldTesla3 = collision3.magField();

        if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
          continue;
        }
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
      }
    } else if (confMixIfTOrVOPartsPresent) {
      Partition<MaskedCollisions> partitionMaskedColT = ((aod::femtodreamcollision::bitmaskTrackOne & maskBit) == maskBit);
      Partition<MaskedCollisions> partitionMaskedColV0 = ((aod::femtodreamcollision::bitmaskTrackThree & maskBit) == maskBit);
      partitionMaskedColT.bindTable(cols);
      partitionMaskedColV0.bindTable(cols);

      for (const auto& [ColWithTrack1, ColWithTrack2, ColWithV0] : soa::combinations(soa::CombinationsBlockFullIndexPolicy(colBinning, confNEventsMix, -1, *partitionMaskedColT.mFiltered, *partitionMaskedColT.mFiltered, *partitionMaskedColV0.mFiltered))) {
        if (ColWithTrack1.globalIndex() == ColWithTrack2.globalIndex() || ColWithTrack1.globalIndex() == ColWithV0.globalIndex() || ColWithTrack2.globalIndex() == ColWithV0.globalIndex()) {
          continue;
        }
        if (ColWithTrack1.globalIndex() > ColWithTrack2.globalIndex()) {
          continue;
        }
        const int multiplicityCol = ColWithTrack1.multNtr();
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({ColWithTrack1.posZ(), multiplicityCol}));

        auto groupPartsOne = selectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithTrack1.globalIndex(), cache);
        auto groupPartsTwo = selectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithTrack2.globalIndex(), cache);
        auto groupPartsThree = selectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithV0.globalIndex(), cache);

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
  PROCESS_SWITCH(FemtoDreamTripletTaskTrackTrackV0PbPb, processMixedEventMasked, "Enable processing mixed events", false);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(const o2::aod::FDCollisions& cols,
                           const soa::Join<o2::aod::FDParticles, const o2::aod::FDMCLabels>& parts,
                           const o2::aod::FDMCParticles&)
  {
    for (const auto& [collision1, collision2, collision3] : soa::selfCombinations(colBinning, confNEventsMix, -1, cols, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      threeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = selectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = selectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsThree = selectedV0sMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

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
  PROCESS_SWITCH(FemtoDreamTripletTaskTrackTrackV0PbPb, processMixedEventMC, "Enable processing mixed events MC", false);

  /// brief process function for to call doMixedEvent with Monte Carlo which has a mask for containing particles or not
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMCMasked(const MaskedCollisions& cols,
                                 const soa::Join<o2::aod::FDParticles, const o2::aod::FDMCLabels>& parts,
                                 const o2::aod::FDMCParticles&)
  {
    if (confMixIfTripletPresent || confMixIfTVOPairPresent) {
      Partition<MaskedCollisions> partitionMaskedCol1 = (confMixIfTripletPresent && (aod::femtodreamcollision::bitmaskTrackTwo & maskBit) == maskBit && (aod::femtodreamcollision::bitmaskTrackThree & maskBit) == maskBit) ||
                                                        (confMixIfTVOPairPresent && (aod::femtodreamcollision::bitmaskTrackOne & maskBit) == maskBit && (aod::femtodreamcollision::bitmaskTrackThree & maskBit) == maskBit);
      partitionMaskedCol1.bindTable(cols);

      for (const auto& [collision1, collision2, collision3] : soa::selfCombinations(colBinning, confNEventsMix, -1, *partitionMaskedCol1.mFiltered, *partitionMaskedCol1.mFiltered, *partitionMaskedCol1.mFiltered)) {
        const int multiplicityCol = collision1.multNtr();
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

        auto groupPartsOne = selectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = selectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
        auto groupPartsThree = selectedV0sMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

        const auto& magFieldTesla1 = collision1.magField();
        const auto& magFieldTesla2 = collision2.magField();
        const auto& magFieldTesla3 = collision3.magField();

        if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
          continue;
        }
        doMixedEvent<true>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
      }
    } else if (confMixIfTOrVOPartsPresent) {
      Partition<MaskedCollisions> partitionMaskedColT = ((aod::femtodreamcollision::bitmaskTrackOne & maskBit) == maskBit);
      Partition<MaskedCollisions> partitionMaskedColV0 = ((aod::femtodreamcollision::bitmaskTrackThree & maskBit) == maskBit);
      partitionMaskedColT.bindTable(cols);
      partitionMaskedColV0.bindTable(cols);

      for (const auto& [ColWithTrack1, ColWithTrack2, ColWithV0] : soa::combinations(soa::CombinationsBlockFullIndexPolicy(colBinning, confNEventsMix, -1, *partitionMaskedColT.mFiltered, *partitionMaskedColT.mFiltered, *partitionMaskedColV0.mFiltered))) {
        if (ColWithTrack1.globalIndex() == ColWithTrack2.globalIndex() || ColWithTrack1.globalIndex() == ColWithV0.globalIndex() || ColWithTrack2.globalIndex() == ColWithV0.globalIndex()) {
          continue;
        }
        if (ColWithTrack1.globalIndex() > ColWithTrack2.globalIndex()) {
          continue;
        }
        const int multiplicityCol = ColWithTrack1.multNtr();
        threeBodyQARegistry.fill(HIST("TripletTaskQA/hMECollisionBins"), colBinning.getBin({ColWithTrack1.posZ(), multiplicityCol}));

        auto groupPartsOne = selectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithTrack1.globalIndex(), cache);
        auto groupPartsTwo = selectedPartsMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithTrack2.globalIndex(), cache);
        auto groupPartsThree = selectedV0sMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, ColWithV0.globalIndex(), cache);

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
  PROCESS_SWITCH(FemtoDreamTripletTaskTrackTrackV0PbPb, processMixedEventMCMasked, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoDreamTripletTaskTrackTrackV0PbPb>(cfgc),
  };
  return workflow;
}
