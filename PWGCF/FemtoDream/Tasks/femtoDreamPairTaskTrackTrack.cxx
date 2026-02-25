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
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamCollisionSelection.h"
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

struct femtoDreamPairTaskTrackTrack {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  /// General options
  struct : ConfigurableGroup {
    std::string prefix = std::string("Option");
    Configurable<bool> IsMC{"IsMC", false, "Enable additional Histogramms in the case of runninger over Monte Carlo"};
    Configurable<bool> Use4D{"Use4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
    Configurable<bool> ExtendedPlots{"ExtendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
    Configurable<float> HighkstarCut{"HighkstarCut", -1., "Set a cut for high k*, above which the pairs are rejected. Set it to -1 to deactivate it"};
    Configurable<bool> SameSpecies{"SameSpecies", false, "Set to true if particle 1 and particle 2 are the same species"};
    Configurable<bool> MixEventWithPairs{"MixEventWithPairs", false, "Only use events that contain particle 1 and partile 2 for the event mixing"};
    Configurable<bool> RandomizePair{"RandomizePair", true, "Randomly mix particle 1 and particle 2 in case both are identical"};
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
    Configurable<float> SphericityMin{"SphericityMin", 0, "Minimum event sphericity"};
  } EventSel;

  Filter EventMultiplicity = aod::femtodreamcollision::multNtr >= EventSel.MultMin && aod::femtodreamcollision::multNtr <= EventSel.MultMax && aod::femtodreamcollision::sphericity >= EventSel.SphericityMin;
  Filter EventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= EventSel.MultPercentileMin && aod::femtodreamcollision::multV0M <= EventSel.MultPercentileMax;

  /// qn&event_plane separator
  FemtoDreamCollisionSelection epCalculator;
  struct : ConfigurableGroup {
    std::string prefix = std::string("EPCal");
    Configurable<bool> do1DFemto{"do1DFemto", false, "Do 1D femtoscopy"};
    Configurable<bool> do3DFemto{"do3DFemto", false, "Do 3D femtoscopy"};
    Configurable<bool> fillFlowQA{"fillFlowQA", false, "Fill QA histos for flow/event-plane related observables"};
    Configurable<bool> storeEvtTrkInfo{"storeEvtTrkInfo", false, "Fill info of track1 and track2 while pariing in divided qn bins"};
    Configurable<bool> doQnSeparation{"doQnSeparation", false, "Do qn separation"};
    Configurable<bool> doEPReClibForMixing{"doEPReClibForMixing", false, "While mixing, using respective event plane for participating particles azimuthal angle caulculation"};
    Configurable<bool> mcQvec{"mcQvec", false, "Enable Q vector table for Monte Carlo"};
    Configurable<std::vector<float>> qnBinSeparator{"qnBinSeparator", std::vector<float>{-999.f, -999.f, -999.f}, "Qn bin separator"};
    Configurable<int> numQnBins{"numQnBins", 10, "Number of qn bins"};
    Configurable<int> qnBinMin{"qnBinMin", 0, "Number of qn bins"};
    Configurable<float> centMax{"centMax", 100.f, "Evt sel: Maximum Centrality cut"};
    Configurable<float> centBinWidth{"centBinWidth", 1.f, "Centrality bin length for qn separator"};
    ConfigurableAxis DKout{"DKout", {500, -2., 2.}, "binning DKout for the 3-D femtoscopy plot: R_out(LCMS) vs mT vs multiplicity percentile vs qnBin vs pait phi wrt EP (set <<do3DFemto>> to true)"};
    ConfigurableAxis DKside{"DKside", {500, -2., 2.}, "binning DKside for the 3-D femtoscopy plot: R_side(LCMS) vs mT vs multiplicity percentile vs qnBin vs pait phi wrt EP (set <<do3DFemto>> to true)"};
    ConfigurableAxis DKlong{"DKlong", {500, -2., 2.}, "binning DKlong for the 3-D femtoscopy plot: R_long(LCMS) vs mT vs multiplicity percentile vs qnBin vs pait phi wrt EP (set <<do3DFemto>> to true)"};
    ConfigurableAxis qnBins{"qnBins", {10, 0, 10}, "binning of qn interval"};
    ConfigurableAxis pairPhiBins{"pairPhiBins", {12, 0., TMath::Pi()}, "binning of pair phi"};
  } EPCal;

  using FilteredCollisions = soa::Filtered<FDCollisions>;
  using FilteredCollision = FilteredCollisions::iterator;
  using FilteredMCCollisions = soa::Filtered<soa::Join<aod::FDCollisions, aod::FDMCCollLabels>>;
  using FilteredMCCollision = FilteredMCCollisions::iterator;
  using FilteredQnCollisions = soa::Filtered<soa::Join<aod::FDCollisions, aod::FDExtQnCollisions, aod::FDExtEPCollisions>>;
  using FilteredQnCollision = FilteredQnCollisions::iterator;
  using FilteredMCQnCollisions = soa::Filtered<soa::Join<aod::FDCollisions, aod::FDMCCollLabels, aod::FDExtQnCollisions, aod::FDExtEPCollisions>>;
  using FilteredMCQnCollision = FilteredMCQnCollisions::iterator;

  using FilteredMaskedCollisions = soa::Filtered<soa::Join<FDCollisions, FDColMasks, FDDownSample>>;
  using FilteredMaskedCollision = FilteredMaskedCollisions::iterator;
  using FilteredMaskedMCCollisions = soa::Filtered<soa::Join<FDCollisions, aod::FDMCCollLabels, FDColMasks, FDDownSample>>;
  using FilteredMaskedMCCollision = FilteredMaskedMCCollisions::iterator;

  femtodreamcollision::BitMaskType BitMask = 0;

  /// Track 1
  struct : ConfigurableGroup {
    std::string prefix = std::string("Track1");
    Configurable<int> PDGCode{"PDGCode", 2212, "PDG code of particle 1 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> CutBit{"CutBit", 3191978, "Selection bit from cutCulator for particle 1 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> TPCBit{"TPCBit", 4, "PID TPC bit from cutCulator for particle 1 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> TPCBit_Reject{"TPCBit_Reject", 0, "PID TPC bit from cutCulator to reject a particle hypothesis for particle 1 (set to 0 to ignore)"};
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
  Partition<aod::FDParticles> PartitionTrk1 =
    (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
    (ncheckbit(aod::femtodreamparticle::cut, Track1.CutBit)) &&
    ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= Track1.PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, Track1.TPCBit) && ((aod::femtodreamparticle::pidcut & Track1.TPCBit_Reject) == 0u), ncheckbit(aod::femtodreamparticle::pidcut, Track1.TPCTOFBit)) &&
    (aod::femtodreamparticle::pt >= Track1.PtMin) &&
    (aod::femtodreamparticle::pt <= Track1.PtMax) &&
    (aod::femtodreamparticle::eta >= Track1.EtaMin) &&
    (aod::femtodreamparticle::eta <= Track1.EtaMax) &&
    ifnode(Option.DCACutPtDep, (nabs(aod::femtodreamparticle::tempFitVar) <= 0.0105f + (0.035f / npow(aod::femtodreamparticle::pt, 1.1f))),
           ((aod::femtodreamparticle::tempFitVar >= Track1.TempFitVarMin) &&
            (aod::femtodreamparticle::tempFitVar <= Track1.TempFitVarMax)));

  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> PartitionMCTrk1 =
    (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
    (ncheckbit(aod::femtodreamparticle::cut, Track1.CutBit)) &&
    ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= Track1.PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, Track1.TPCBit) && ((aod::femtodreamparticle::pidcut & Track1.TPCBit_Reject) == 0u), ncheckbit(aod::femtodreamparticle::pidcut, Track1.TPCTOFBit)) &&
    (aod::femtodreamparticle::pt >= Track1.PtMin) &&
    (aod::femtodreamparticle::pt <= Track1.PtMax) &&
    (aod::femtodreamparticle::eta >= Track1.EtaMin) &&
    (aod::femtodreamparticle::eta < Track1.EtaMax) &&
    ifnode(Option.DCACutPtDep, (nabs(aod::femtodreamparticle::tempFitVar) <= 0.0105f + (0.035f / npow(aod::femtodreamparticle::pt, 1.1f))),
           ((aod::femtodreamparticle::tempFitVar >= Track1.TempFitVarMin) &&
            (aod::femtodreamparticle::tempFitVar <= Track1.TempFitVarMax)));

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Track 2
  struct : ConfigurableGroup {
    std::string prefix = std::string("Track2");
    Configurable<int> PDGCode{"PDGCode", 2212, "PDG code of particle 2 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> CutBit{"CutBit", 3191978, "Selection bit from cutCulator for particle 2 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> TPCBit{"TPCBit", 4, "PID TPC bit from cutCulator for particle 2 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> TPCBit_Reject{"TPCBit_Reject", 0, "PID TPC bit from cutCulator to reject a particle hypothesis for particle 2 (set to 0 to ignore)"};
    Configurable<femtodreamparticle::cutContainerType> TPCTOFBit{"TPCTOFBit", 2, "PID TPCTOF bit from cutCulator for particle 2 (Track)"};
    Configurable<float> PIDThres{"PIDThres", 0.75, "Momentum threshold for PID selection for particle 2 (Track)"};
    Configurable<float> PtMin{"PtMin", 0., "Minimum pT of particle 2 (Track)"};
    Configurable<float> PtMax{"PtMax", 999., "Maximum pT of particle 2 (Track)"};
    Configurable<float> EtaMin{"EtaMin", -10., "Minimum eta of particle 2 (Track)"};
    Configurable<float> EtaMax{"EtaMax", 10., "Maximum eta of particle 2 (Track)"};
    Configurable<float> TempFitVarMin{"TempFitVarMin", -10., "Minimum DCAxy of partricle 1 (Track)"};
    Configurable<float> TempFitVarMax{"TempFitVarMax", 10., "Maximum DCAxy of partricle 1 (Track)"};
  } Track2;

  /// Partition for track 2
  Partition<aod::FDParticles> PartitionTrk2 =
    (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
    (ncheckbit(aod::femtodreamparticle::cut, Track2.CutBit)) &&
    ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= Track2.PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, Track2.TPCBit) && ((aod::femtodreamparticle::pidcut & Track2.TPCBit_Reject) == 0u), ncheckbit(aod::femtodreamparticle::pidcut, Track2.TPCTOFBit)) &&
    (aod::femtodreamparticle::pt >= Track2.PtMin) &&
    (aod::femtodreamparticle::pt <= Track2.PtMax) &&
    (aod::femtodreamparticle::eta >= Track2.EtaMin) &&
    (aod::femtodreamparticle::eta <= Track2.EtaMax) &&
    ifnode(Option.DCACutPtDep, (nabs(aod::femtodreamparticle::tempFitVar) <= 0.0105f + (0.035f / npow(aod::femtodreamparticle::pt, 1.1f))),
           ((aod::femtodreamparticle::tempFitVar >= Track2.TempFitVarMin) &&
            (aod::femtodreamparticle::tempFitVar <= Track2.TempFitVarMax)));

  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> PartitionMCTrk2 =
    (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
    (ncheckbit(aod::femtodreamparticle::cut, Track2.CutBit)) &&
    ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= Track2.PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, Track2.TPCBit) && ((aod::femtodreamparticle::pidcut & Track2.TPCBit_Reject) == 0u), ncheckbit(aod::femtodreamparticle::pidcut, Track2.TPCTOFBit)) &&
    (aod::femtodreamparticle::pt >= Track2.PtMin) &&
    (aod::femtodreamparticle::pt <= Track2.PtMax) &&
    (aod::femtodreamparticle::eta > Track2.EtaMin) &&
    (aod::femtodreamparticle::eta < Track2.EtaMax) &&
    ifnode(Option.DCACutPtDep, (nabs(aod::femtodreamparticle::tempFitVar) <= 0.0105f + (0.035f / npow(aod::femtodreamparticle::pt, 1.1f))),
           ((aod::femtodreamparticle::tempFitVar >= Track2.TempFitVarMin) &&
            (aod::femtodreamparticle::tempFitVar <= Track2.TempFitVarMax)));

  /// Histogramming for track 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// Binning configurables
  struct : ConfigurableGroup {
    std::string prefix = "Binning";
    ConfigurableAxis TempFitVar{"TempFitVar", {300, -0.15, 0.15}, "Binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis TrackpT{"TrackpT", {20, 0.5, 4.05}, "pT binning for pT vs. TempFitVar plot"};
    ConfigurableAxis pT{"pT", {20, 0.5, 4.05}, "pT binning for extended plots"};
    ConfigurableAxis kstar{"kstar", {1500, 0., 6.}, "kstar binning"};
    ConfigurableAxis kT{"kT", {150, 0., 9.}, "kT binning"};
    ConfigurableAxis mT{"mT", {225, 0., 7.5}, "mT binning"};
    ConfigurableAxis multTempFit{"multTempFit", {1, 0, 1}, "Multiplicity Binning for the TempFitVar plot"};
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
    ConfigurableAxis QnMixBins{"QnMixBins", {VARIABLE_WIDTH, 0.50f, 68.50f, 100.50f, 126.50f, 151.50f, 176.50f, 203.50f, 232.50f, 269.50f, 322.50f, 833.50f}, "Mixing bins - qn-value"};
    ConfigurableAxis EPMixBins{"EPMixBins", {VARIABLE_WIDTH, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180.}, "Mixing bins - event plane (deg)"};
    Configurable<int> Depth{"Depth", 5, "Number of events for mixing"};
    Configurable<int> Policy{"Policy", 0, "Binning policy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both, 3: multipliciy percentile and qn value, 4: multipliciy percentile and event plane"};
  } Mixing;

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{Mixing.VztxMixBins, Mixing.MultMixBins}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{Mixing.VztxMixBins, Mixing.MultPercentileMixBins}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{Mixing.VztxMixBins, Mixing.MultMixBins, Mixing.MultPercentileMixBins}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M, aod::femtodreamcollision::QnVal> colBinningMultPercentileqn{{Mixing.VztxMixBins, Mixing.MultPercentileMixBins, Mixing.QnMixBins}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M, aod::femtodreamcollision::EventPlane> colBinningMultPercentileEP{{Mixing.VztxMixBins, Mixing.MultPercentileMixBins, Mixing.EPMixBins}, true};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejectionSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejectionME;

  // Container for correlation functions in devided qn bins
  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventQnCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventQnCont;

  /// Histogram output
  HistogramRegistry Registry{"Output", {}, OutputObjHandlingPolicy::AnalysisObject};

  TRandom3* random;

  void init(InitContext& context)
  {

    // setup columnpolicy for binning
    colBinningMult = {{Mixing.VztxMixBins, Mixing.MultMixBins}, true};
    colBinningMultPercentile = {{Mixing.VztxMixBins, Mixing.MultPercentileMixBins}, true};
    colBinningMultMultPercentile = {{Mixing.VztxMixBins, Mixing.MultMixBins, Mixing.MultPercentileMixBins}, true};
    colBinningMultPercentileqn = {{Mixing.VztxMixBins, Mixing.MultPercentileMixBins, Mixing.QnMixBins}, true};
    colBinningMultPercentileEP = {{Mixing.VztxMixBins, Mixing.MultPercentileMixBins, Mixing.EPMixBins}, true};

    if (Option.RandomizePair.value) {
      random = new TRandom3(0);
    }
    eventHisto.init(&Registry, Option.IsMC);
    trackHistoPartOne.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.TrackpT, Option.Dummy, Option.Dummy, Binning.TempFitVar, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.IsMC, Track1.PDGCode);
    if (!Option.SameSpecies) {
      trackHistoPartTwo.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.TrackpT, Option.Dummy, Option.Dummy, Binning.TempFitVar, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.IsMC, Track2.PDGCode);
    }

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
    sameEventCont.setPDGCodes(Track1.PDGCode, Track2.PDGCode);
    mixedEventCont.setPDGCodes(Track1.PDGCode, Track2.PDGCode);
    pairCleaner.init(&Registry);
    if (Option.CPROn.value) {
      pairCloseRejectionSE.init(&Registry, &Registry, Option.CPRdeltaPhiMax.value, Option.CPRdeltaEtaMax.value, Option.CPRPlotPerRadii.value, 1, Option.CPROld.value);
      pairCloseRejectionME.init(&Registry, &Registry, Option.CPRdeltaPhiMax.value, Option.CPRdeltaEtaMax.value, Option.CPRPlotPerRadii.value, 2, Option.CPROld.value);
    }

    if (EPCal.do1DFemto) {
      sameEventQnCont.init_EP(&Registry,
                              Binning4D.kstar, Binning4D.mT, Binning4D.multPercentile, EPCal.doQnSeparation ? EPCal.qnBins : EPCal.pairPhiBins, EPCal.doQnSeparation ? "qnBin" : "#phi_{pair}-#Psi_{EP} (rad)",
                              Option.IsMC, Option.HighkstarCut);
      mixedEventQnCont.init_EP(&Registry,
                               Binning4D.kstar, Binning4D.mT, Binning4D.multPercentile, EPCal.doQnSeparation ? EPCal.qnBins : EPCal.pairPhiBins, EPCal.doQnSeparation ? "qnBin" : "#phi_{pair}-#Psi_{EP} (rad)",
                               Option.IsMC, Option.HighkstarCut);
      sameEventQnCont.setPDGCodes(Track1.PDGCode, Track2.PDGCode);
      mixedEventQnCont.setPDGCodes(Track1.PDGCode, Track2.PDGCode);
      if (EPCal.fillFlowQA) {
        epCalculator.initEPQA(&Registry);
      }
    }

    if (EPCal.do3DFemto) {
      sameEventQnCont.init_3Dqn(&Registry, EPCal.DKout, EPCal.DKside, EPCal.DKlong,
                                Binning4D.mT, Binning4D.multPercentile, Option.IsMC, EPCal.qnBins, EPCal.pairPhiBins, Option.SmearingByOrigin);
      mixedEventQnCont.init_3Dqn(&Registry, EPCal.DKout, EPCal.DKside, EPCal.DKlong,
                                 Binning4D.mT, Binning4D.multPercentile, Option.IsMC, EPCal.qnBins, EPCal.pairPhiBins, Option.SmearingByOrigin);
      sameEventQnCont.setPDGCodes(Track1.PDGCode, Track2.PDGCode);
      mixedEventQnCont.setPDGCodes(Track1.PDGCode, Track2.PDGCode);
      if (EPCal.fillFlowQA) {
        epCalculator.initEPQA(&Registry);
      }
    }

    // get bit for the collision mask
    std::bitset<8 * sizeof(femtodreamcollision::BitMaskType)> mask;
    int index = 0;
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.find("femto-dream-pair-task-track-track") != std::string::npos) {
        if (containsNameValuePair(device.options, "Option.DCACutPtDep", Option.DCACutPtDep.value) &&
            containsNameValuePair(device.options, "Option.SameSpecies", Option.SameSpecies.value) &&
            containsNameValuePair(device.options, "Track1.CutBit", Track1.CutBit.value) &&
            containsNameValuePair(device.options, "Track1.TPCBit", Track1.TPCBit.value) &&
            containsNameValuePair(device.options, "Track1.TPCTOFBit", Track1.TPCTOFBit.value) &&
            containsNameValuePair(device.options, "Track1.PIDThres", Track1.PIDThres.value) &&
            containsNameValuePair(device.options, "Track1.PtMin", Track1.PtMin.value) &&
            containsNameValuePair(device.options, "Track1.PtMax", Track1.PtMax.value) &&
            containsNameValuePair(device.options, "Track1.EtaMin", Track1.EtaMin.value) &&
            containsNameValuePair(device.options, "Track1.EtaMax", Track1.EtaMax.value) &&
            containsNameValuePair(device.options, "Track1.TempFitVarMin", Track1.TempFitVarMin.value) &&
            containsNameValuePair(device.options, "Track1.TempFitVarMax", Track1.TempFitVarMax.value) &&
            containsNameValuePair(device.options, "Track2.CutBit", Track2.CutBit.value) &&
            containsNameValuePair(device.options, "Track2.TPCBit", Track2.TPCBit.value) &&
            containsNameValuePair(device.options, "Track2.TPCTOFBit", Track2.TPCTOFBit.value) &&
            containsNameValuePair(device.options, "Track2.PIDThres", Track2.PIDThres.value) &&
            containsNameValuePair(device.options, "Track2.PtMin", Track2.PtMin.value) &&
            containsNameValuePair(device.options, "Track2.PtMax", Track2.PtMax.value) &&
            containsNameValuePair(device.options, "Track2.EtaMin", Track2.EtaMin.value) &&
            containsNameValuePair(device.options, "Track2.EtaMax", Track2.EtaMax.value) &&
            containsNameValuePair(device.options, "Track2.TempFitVarMin", Track2.TempFitVarMin.value) &&
            containsNameValuePair(device.options, "Track2.TempFitVarMax", Track2.TempFitVarMax.value)) {
          mask.set(index);
          BitMask = static_cast<femtodreamcollision::BitMaskType>(mask.to_ulong());
          LOG(info) << "Configuration matched for device: " << device.name;
          LOG(info) << "Bitmask for collisions: " << mask.to_string();
          break;
        } else {
          index++;
        }
      }
    }
    if ((doprocessSameEvent && doprocessSameEventMasked) ||
        (doprocessSameEvent && doprocessSameEventEP) ||
        (doprocessSameEventMasked && doprocessSameEventEP) ||
        (doprocessMixedEvent && doprocessMixedEventMasked) ||
        (doprocessMixedEvent && doprocessMixedEventEP) ||
        (doprocessMixedEventMasked && doprocessMixedEventEP) ||
        (doprocessSameEventMC && doprocessSameEventMCMasked) ||
        (doprocessSameEventMC && doprocessSameEventEPMC) ||
        (doprocessMixedEventMC && doprocessMixedEventMCMasked) ||
        (doprocessMixedEventMC && doprocessMixedEventEPMC) ||
        (doprocessMixedEventMCMasked && doprocessMixedEventEPMC)) {
      LOG(fatal) << "Normal and masked processing cannot be activated simultaneously!";
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
  void doSameEvent(PartitionType SliceTrk1, PartitionType SliceTrk2, PartType parts, Collision col)
  {
    for (auto& part : SliceTrk1) {
      trackHistoPartOne.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
    }

    if (!Option.SameSpecies.value) {
      for (auto& part : SliceTrk2) {
        trackHistoPartTwo.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }

    /// Now build the combinations
    float rand = 0.;
    if (Option.SameSpecies.value) {
      for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(SliceTrk1, SliceTrk2))) {
        if (Option.CPROn.value) {
          if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (Option.RandomizePair.value) {
          rand = random->Rndm();
        }
        if (rand <= 0.5) {
          sameEventCont.setPair<isMC>(p1, p2, col.multNtr(), col.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.SmearingByOrigin);
        } else {
          sameEventCont.setPair<isMC>(p2, p1, col.multNtr(), col.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.SmearingByOrigin);
        }
      }
    } else {
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceTrk2))) {
        if (Option.CPROn.value) {
          if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
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
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceTrk2 = PartitionTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceTrk1.size() == 0 && SliceTrk2.size() == 0) {
      return;
    }
    doSameEvent<false>(SliceTrk1, SliceTrk2, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEvent, "Enable processing same event", true);

  void processSameEventMasked(FilteredMaskedCollision& col, o2::aod::FDParticles& parts)
  {
    if (Option.SameSpecies.value) {
      if ((col.bitmaskTrackOne() & BitMask) != BitMask) {
        return;
      }
    } else {
      if ((col.bitmaskTrackOne() & BitMask) != BitMask || (col.bitmaskTrackTwo() & BitMask) != BitMask) {
        return;
      }
    }
    fillCollision<false>(col);
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceTrk2 = PartitionTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<false>(SliceTrk1, SliceTrk2, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEventMasked, "Enable processing same event with masks", false);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(FilteredMCCollision& col,
                          o2::aod::FDMCCollisions&,
                          soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                          o2::aod::FDMCParticles&)
  {
    fillCollision<true>(col);
    auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceMCTrk2 = PartitionMCTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceMCTrk1.size() == 0 && SliceMCTrk2.size() == 0) {
      return;
    }
    doSameEvent<true>(SliceMCTrk1, SliceMCTrk2, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  void processSameEventMCMasked(FilteredMaskedMCCollision& col, o2::aod::FDMCCollisions&, soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                                o2::aod::FDMCParticles&)
  {
    if (Option.SameSpecies.value) {
      if ((col.bitmaskTrackOne() & BitMask) != BitMask) {
        return;
      }
    } else {
      if ((col.bitmaskTrackOne() & BitMask) != BitMask && (col.bitmaskTrackTwo() & BitMask) != BitMask) {
        return;
      }
    }
    fillCollision<true>(col);
    auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceMCTrk2 = PartitionMCTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<true>(SliceMCTrk1, SliceMCTrk2, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEventMCMasked, "Enable processing same event for Monte Carlo with masked collisions", false);

  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent_NotMasked(CollisionType& cols, PartType& parts, PartitionType& part1, PartitionType& part2, BinningType policy)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(policy, Mixing.Depth.value, -1, cols, cols)) {
      auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceTrk2 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      if (SliceTrk1.size() == 0 || SliceTrk2.size() == 0) {
        continue;
      }
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceTrk2))) {
        if (Option.CPROn.value) {
          if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
            continue;
          }
        }
        mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.SmearingByOrigin);
      }
    }
  }

  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent_Masked(CollisionType& cols, PartType& parts, PartitionType& part1, PartitionType& part2, BinningType policy)
  {
    if (!Option.SameSpecies.value && !Option.MixEventWithPairs.value) {
      // If the two particles are not the same species and the events which are mixed should contain at least one particle of interest, create two paritition of collisions that contain at least one of the two particle of interest and mix them
      // Make sure there is a check that we do not mix a event with itself in case it contains both partilces
      Partition<CollisionType> PartitionMaskedCol1 = (aod::femtodreamcollision::bitmaskTrackOne & BitMask) == BitMask && aod::femtodreamcollision::downsample == true;
      PartitionMaskedCol1.bindTable(cols);
      Partition<CollisionType> PartitionMaskedCol2 = (aod::femtodreamcollision::bitmaskTrackTwo & BitMask) == BitMask && aod::femtodreamcollision::downsample == true;
      PartitionMaskedCol2.bindTable(cols);
      // use *Partition.mFiltered when passing the partition to mixing object
      // there is an issue when the partition is passed directly
      // workaround for now, change back once it is fixed
      for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockUpperIndexPolicy(policy, Mixing.Depth.value, -1, *PartitionMaskedCol1.mFiltered, *PartitionMaskedCol2.mFiltered))) {
        // make sure that tracks in the same events are not mixed
        if (collision1.globalIndex() == collision2.globalIndex()) {
          continue;
        }
        auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto SliceTrk2 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);

        for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceTrk2))) {
          if (Option.CPROn.value) {
            if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
              continue;
            }
          }
          mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.SmearingByOrigin);
        }
      }
    } else {
      // In the other case where the two particles are not the same species and we do not mix event with pairs,  we only need to define one partition of collisions and make self combinations

      // define a lambda function for the mixing with selfCombinations policy
      auto MixEvents = [policy, &part1, &part2, parts, this](auto& partition) {
        for (auto const& [collision1, collision2] : selfCombinations(policy, Mixing.Depth.value, -1, *partition.mFiltered, *partition.mFiltered)) {
          auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
          auto SliceTrk2 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
          for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceTrk2))) {
            if (Option.CPROn.value) {
              if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
                continue;
              }
            }
            mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), Option.Use4D.value, Option.ExtendedPlots.value, Option.SmearingByOrigin.value);
          }
        }
      };
      if (Option.SameSpecies.value && Option.MixEventWithPairs.value) {
        // in case of mixing the same species and events that contain pairs, check for bitmask of particle two
        // when same species is set to true the bit of particle two is only set if the event contains at least two selected particles
        Partition<CollisionType> PartitionMaskedCol1 = ncheckbit(aod::femtodreamcollision::bitmaskTrackTwo, BitMask) && aod::femtodreamcollision::downsample == true;
        PartitionMaskedCol1.bindTable(cols);
        MixEvents(PartitionMaskedCol1);
      } else if (Option.SameSpecies.value && !Option.MixEventWithPairs.value) {
        // in case of mixing the same species and events that contain at least one selected paritcle, check for bitmask of particle one
        Partition<CollisionType> PartitionMaskedCol1 = ncheckbit(aod::femtodreamcollision::bitmaskTrackOne, BitMask) && aod::femtodreamcollision::downsample == true;
        PartitionMaskedCol1.bindTable(cols);
        MixEvents(PartitionMaskedCol1);
      } else if (!Option.SameSpecies.value && Option.MixEventWithPairs.value) {
        // in case of mixing different species and events that contain a pair of selected paritcles, check for both bitmasks of paritcle one and particle two
        Partition<CollisionType> PartitionMaskedCol1 = ncheckbit(aod::femtodreamcollision::bitmaskTrackOne, BitMask) && ncheckbit(aod::femtodreamcollision::bitmaskTrackTwo, BitMask) && aod::femtodreamcollision::downsample == true;
        PartitionMaskedCol1.bindTable(cols);
        MixEvents(PartitionMaskedCol1);
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
        doMixedEvent_NotMasked<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_NotMasked<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_NotMasked<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEvent, "Enable processing mixed events", true);

  void processMixedEventMasked(FilteredMaskedCollisions& cols, o2::aod::FDParticles& parts)
  {
    switch (Mixing.Policy.value) {
      case static_cast<int>(femtodreamcollision::kMult):
        doMixedEvent_Masked<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_Masked<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_Masked<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEventMasked, "Enable processing mixed events", false);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(FilteredMCCollisions& cols, o2::aod::FDMCCollisions&, soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts, o2::aod::FDMCParticles&)
  {
    switch (Mixing.Policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent_NotMasked<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_NotMasked<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_NotMasked<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEventMC, "Enable processing mixed events MC", false);

  void processMixedEventMCMasked(FilteredMaskedMCCollisions& cols, o2::aod::FDMCCollisions&, soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts, o2::aod::FDMCParticles&)
  {
    switch (Mixing.Policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent_Masked<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_Masked<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_Masked<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEventMCMasked, "Enable processing mixed events MC with masked collisions", false);

  /// This function processes the same event in divided qn bins
  /// col.multV0M() get the event centrality from ft0c for PbPb data
  template <bool isMC, typename PartitionType, typename PartType, typename Collision>
  void doSameEventEP(PartitionType SliceTrk1, PartitionType SliceTrk2, PartType parts, Collision col)
  {
    if (EPCal.storeEvtTrkInfo) {
      for (auto& part : SliceTrk1) {
        trackHistoPartOne.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }

      if (!Option.SameSpecies.value) {
        for (auto& part : SliceTrk2) {
          trackHistoPartTwo.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        }
      }
    }

    float myEP = -999.f;
    int myqnBin = -999;

    if (!isMC || EPCal.mcQvec){
      myEP = TMath::DegToRad() * col.eventPlane();
      if (EPCal.doQnSeparation || EPCal.do3DFemto) {
        myqnBin = epCalculator.myqnBin(col.multV0M(), EPCal.centMax, EPCal.fillFlowQA, EPCal.qnBinSeparator, col.qnVal(), EPCal.numQnBins, EPCal.centBinWidth);
        if (myqnBin < EPCal.qnBinMin || myqnBin > EPCal.numQnBins) {
          myqnBin = -999;
        }
      }
    } else {
      myEP = 0.f;
      myqnBin = 0;
    }

    /// Now build the combinations
    float rand = 0.;
    if (Option.SameSpecies.value) {
      for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(SliceTrk1, SliceTrk2))) {
        if (Option.CPROn.value) {
          if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (Option.RandomizePair.value) {
          rand = random->Rndm();
        }
        if (rand <= 0.5) {
          if (EPCal.do1DFemto) {
            sameEventQnCont.setPair_EP<isMC>(p1, p2, col.multV0M(), EPCal.doQnSeparation, EPCal.doQnSeparation ? myqnBin + 0.f : myEP);
          }
          if (EPCal.do3DFemto) {
            sameEventQnCont.setPair_3Dqn<isMC>(p1, p2, col.multV0M(), Option.SameSpecies.value, myqnBin + 0.f, myEP, Option.SmearingByOrigin);
          } 
        } else {
          if (EPCal.do1DFemto) {
            sameEventQnCont.setPair_EP<isMC>(p2, p1, col.multV0M(), EPCal.doQnSeparation, EPCal.doQnSeparation ? myqnBin + 0.f : myEP);
          }
          if (EPCal.do3DFemto) {
            sameEventQnCont.setPair_3Dqn<isMC>(p2, p1, col.multV0M(), Option.SameSpecies.value, myqnBin + 0.f, myEP, Option.SmearingByOrigin);
          }
        }
      }
    } else {
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceTrk2))) {
        if (Option.CPROn.value) {
          if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (EPCal.do1DFemto) {
          sameEventQnCont.setPair_EP<isMC>(p1, p2, col.multV0M(), EPCal.doQnSeparation, EPCal.doQnSeparation ? myqnBin + 0.f : myEP);
        }
        if (EPCal.do3DFemto) {
          sameEventQnCont.setPair_3Dqn<isMC>(p1, p2, col.multV0M(), Option.SameSpecies.value, myEP, myqnBin, Option.SmearingByOrigin);
        }
      }
    }
  }

  /// process function for to call doSameEventEP with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEventEP(FilteredQnCollision& col, o2::aod::FDParticles& parts)
  {
    if (EPCal.storeEvtTrkInfo) {
      fillCollision<false>(col);
    }
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceTrk2 = PartitionTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceTrk1.size() == 0 && SliceTrk2.size() == 0) {
      return;
    }
    if (EPCal.do1DFemto || EPCal.do3DFemto) {
      doSameEventEP<false>(SliceTrk1, SliceTrk2, parts, col);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEventEP, "Enable processing same event wrt azimuthal angle and event-plane ", false);

  /// process function for to call doSameEventEP with MC Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEventEPMC(FilteredMCQnCollision& col, 
                            o2::aod::FDMCCollisions&,
                            soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                            o2::aod::FDMCParticles&)
  {
    if (EPCal.storeEvtTrkInfo) {
      fillCollision<true>(col);
    }    
    auto SliceTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceTrk2 = PartitionMCTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceTrk1.size() == 0 && SliceTrk2.size() == 0) {
      return;
    }
    if (EPCal.do1DFemto || EPCal.do3DFemto) {
      doSameEventEP<true>(SliceTrk1, SliceTrk2, parts, col);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEventEPMC, "Enable processing same event of 3D for Monte Carlo", false);

  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent_NotMaskedEP(CollisionType& cols, PartType& parts, PartitionType& part1, PartitionType& part2, BinningType policy)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(policy, Mixing.Depth.value, -1, cols, cols)) {
      auto SliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceTrk2 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      if (SliceTrk1.size() == 0 || SliceTrk2.size() == 0) {
        continue;
      }

      auto myEP_event1 = -999.f;
      auto myEP_event2 = -999.f;

      if (!isMC || EPCal.mcQvec){
        myEP_event1 = TMath::DegToRad() * collision1.eventPlane();
        myEP_event2 = TMath::DegToRad() * collision2.eventPlane();
      } 
      else{
        myEP_event1 = 0.f;
        myEP_event2 = 0.f;
      }


      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceTrk2))) {
        if (Option.CPROn.value) {
          if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
            continue;
          }
        }
        if (EPCal.doEPReClibForMixing) {
          if (EPCal.do1DFemto) {
            if (EPCal.doQnSeparation)
              mixedEventQnCont.setPair_EP<isMC>(p1, p2, collision1.multV0M(), EPCal.doQnSeparation, 0.f);
            else
              mixedEventQnCont.setPair_EP<isMC>(p1, p2, collision1.multV0M(), EPCal.doQnSeparation, myEP_event1, myEP_event2);
          }
          if (EPCal.do3DFemto) {
            mixedEventQnCont.setPair_3Dqn<isMC>(p1, p2, collision1.multV0M(), Option.SameSpecies.value, 0.f, myEP_event1, myEP_event2, Option.SmearingByOrigin);
          }
        } else {
          if (EPCal.do1DFemto)
            mixedEventQnCont.setPair_EP<isMC>(p1, p2, collision1.multV0M(), EPCal.doQnSeparation, EPCal.doQnSeparation ? 0.f : myEP_event1);
          if (EPCal.do3DFemto) {
            mixedEventQnCont.setPair_3Dqn<isMC>(p1, p2, collision1.multV0M(), Option.SameSpecies.value, 0.f, myEP_event1, Option.SmearingByOrigin);
          }
        }
      }
    }
  };

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoDreamParticleTable
  void processMixedEventEP(FilteredQnCollisions& cols, o2::aod::FDParticles& parts)
  {
    switch (Mixing.Policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent_NotMaskedEP<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_NotMaskedEP<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_NotMaskedEP<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMultMultPercentile);
        break;
      case femtodreamcollision::kMultPercentileQn:
        doMixedEvent_NotMaskedEP<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMultPercentileqn);
        break;
      case femtodreamcollision::kMultPercentileEP:
        doMixedEvent_NotMaskedEP<false>(cols, parts, PartitionTrk1, PartitionTrk2, colBinningMultPercentileEP);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEventEP, "Enable processing mixed events wrt azimuthal angle and event-plane", false);

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoDreamParticleTable
  void processMixedEventEPMC(FilteredMCQnCollisions& cols, o2::aod::FDMCCollisions&, soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts, o2::aod::FDMCParticles&)
  {
    switch (Mixing.Policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent_NotMaskedEP<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent_NotMaskedEP<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent_NotMaskedEP<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMultMultPercentile);
        break;
      case femtodreamcollision::kMultPercentileQn:
        doMixedEvent_NotMaskedEP<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMultPercentileqn);
        break;
      case femtodreamcollision::kMultPercentileEP:
        doMixedEvent_NotMaskedEP<true>(cols, parts, PartitionMCTrk1, PartitionMCTrk2, colBinningMultPercentileEP);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEventEPMC, "Enable processing mixed events wrt azimuthal angle and event-plane for Monte Carlo", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskTrackTrack>(cfgc),
  };
  return workflow;
}
