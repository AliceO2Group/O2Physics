// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file baryonPairInJets.cxx
/// \brief Identified proton-pair correlations in charged jets
/// \author Hanseo Park <hanseo.park@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;

using BaryonTracks = soa::Join<aod::Tracks,
                               aod::TracksExtra,
                               aod::TracksDCA,
                               aod::TrackSelection,
                               aod::pidTPCFullPi,
                               aod::pidTOFFullPi,
                               aod::pidTPCFullKa,
                               aod::pidTOFFullKa,
                               aod::pidTPCFullPr,
                               aod::pidTOFFullPr,
                               aod::pidTOFbeta>;
using BaryonTracksMC = soa::Join<BaryonTracks, aod::McTrackLabels>;
using JetTracksWithOriginal = soa::Join<aod::JetTracks, aod::JTrackPIs>;
using JetTracksWithOriginalAndMC = soa::Join<aod::JetTracks, aod::JTrackPIs, aod::JMcTrackLbs>;
using ChargedJetsWithConstituents = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;
using ChargedMCDJetsWithConstituents = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;
using ChargedMCPJetsWithConstituents = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;
using MatchedChargedMCDJets = soa::Join<ChargedMCDJetsWithConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>;
using MatchedChargedMCPJets = soa::Join<ChargedMCPJetsWithConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;

struct BaryonPairInJets {
  enum class PairType : uint8_t {
    SameSign = 0,
    OppositeSign = 1
  };

  enum class ChargePairType : uint8_t {
    PositivePositive = 0,
    NegativeNegative = 1,
    OppositeSign = 2
  };

  enum class MixedReference : uint8_t {
    Nominal = 0,
    Loose = 1,
    Random = 2
  };

  enum class ProtonPidMethod : int {
    TpcOnly = 0,
    TpcTofCircular = 1,
    TpcTofRectangular = 2
  };

  enum class ControlSpecies : uint8_t {
    Pion = 0,
    Kaon = 1
  };

  enum class EfficiencyStage : uint8_t {
    Generated = 0,
    TrackSelected = 1,
    PidSelected = 2
  };

  // Species-independent selected-particle representation. Keep the full
  // kinematics needed by pair and future jet-fragmentation observables.
  struct SelectedParticle {
    float pt;
    float eta;
    float phi;
    int8_t sign;
  };

  struct SelectedParticles {
    std::vector<SelectedParticle> baryons;
    std::vector<SelectedParticle> pions;
    std::vector<SelectedParticle> kaons;
  };

  struct : ConfigurableGroup {
    Configurable<float> vertexZMax{"vertexZMax", 10.f, "Maximum absolute collision z vertex (cm)"};
    Configurable<std::string> eventSelections{"eventSelections", "sel8+NoITSROFrameBorder+NoTimeFrameBorder+NoSameBunchPileup+IsGoodZvtxFT0vsPV", "Event selections"};
    Configurable<bool> applyRctSelection{"applyRctSelection", true, "Apply the PWGJE hadron-PID RCT selection"};
    Configurable<std::string> rctLabel{"rctLabel", "CBT_hadronPID", "RCT selection label"};
    Configurable<int> occupancyMin{"occupancyMin", -999999, "Minimum track occupancy in time range"};
    Configurable<int> occupancyMax{"occupancyMax", 999999, "Maximum track occupancy in time range"};
  } eventCuts;

  struct : ConfigurableGroup {
    Configurable<float> jetRadius{"jetRadius", 0.4f, "Charged-jet resolution parameter"};
    Configurable<float> jetPtMin{"jetPtMin", 5.f, "Minimum charged-jet transverse momentum (GeV/c); broad input range"};
    Configurable<float> jetPtMax{"jetPtMax", 140.f, "Maximum charged-jet transverse momentum (GeV/c); broad input range"};
    Configurable<float> jetEtaMax{"jetEtaMax", 0.5f, "Maximum absolute charged-jet pseudorapidity"};
    Configurable<float> areaFractionMin{"areaFractionMin", -1.f, "Minimum jet area divided by pi R squared; negative disables the cut"};
  } jetCuts;

  struct : ConfigurableGroup {
    Configurable<float> protonPtMin{"protonPtMin", 0.5f, "Minimum proton transverse momentum (GeV/c)"};
    Configurable<float> protonPtMax{"protonPtMax", 2.5f, "Maximum proton transverse momentum (GeV/c)"};
    Configurable<float> protonEtaMax{"protonEtaMax", 0.8f, "Maximum absolute proton pseudorapidity"};
    Configurable<bool> requirePvContributor{"requirePvContributor", false, "Require the track to be a PV contributor"};
    Configurable<int> tpcCrossedRowsMin{"tpcCrossedRowsMin", 70, "Minimum number of crossed TPC rows"};
    Configurable<float> tpcChi2Max{"tpcChi2Max", 4.f, "Maximum TPC chi2 per cluster"};
    Configurable<float> itsChi2Max{"itsChi2Max", 36.f, "Maximum ITS chi2 per cluster"};
    Configurable<float> dcaXYMax{"dcaXYMax", 0.2f, "Maximum absolute transverse DCA (cm)"};
    Configurable<float> dcaZMax{"dcaZMax", 0.1f, "Maximum absolute longitudinal DCA (cm)"};
  } trackCuts;

  struct : ConfigurableGroup {
    Configurable<int> method{"method", static_cast<int>(ProtonPidMethod::TpcTofCircular), "PID method: 0 TPC only, 1 circular TPC-TOF, 2 rectangular TPC-TOF"};
    Configurable<float> tofPtThreshold{"tofPtThreshold", 0.75f, "Require TOF above this transverse momentum (GeV/c)"};
    Configurable<float> tpcNSigmaMax{"tpcNSigmaMax", 3.f, "Maximum absolute proton TPC n-sigma"};
    Configurable<float> tofNSigmaMax{"tofNSigmaMax", 3.f, "Maximum absolute proton TOF n-sigma"};
    Configurable<float> combinedNSigmaMax{"combinedNSigmaMax", 3.f, "Maximum circular TPC-TOF proton n-sigma"};
    Configurable<bool> rejectPionKaon{"rejectPionKaon", false, "Reject candidates compatible with the pion or kaon hypothesis"};
    Configurable<float> rejectionNSigma{"rejectionNSigma", 3.f, "Pion and kaon rejection threshold"};
    Configurable<float> pionTpcNSigmaMax{"pionTpcNSigmaMax", 3.f, "Maximum absolute pion TPC n-sigma"};
    Configurable<float> pionTofNSigmaMax{"pionTofNSigmaMax", 3.f, "Maximum absolute pion TOF n-sigma"};
    Configurable<float> pionCombinedNSigmaMax{"pionCombinedNSigmaMax", 3.f, "Maximum circular TPC-TOF pion n-sigma"};
    Configurable<float> kaonTpcNSigmaMax{"kaonTpcNSigmaMax", 3.f, "Maximum absolute kaon TPC n-sigma"};
    Configurable<float> kaonTofNSigmaMax{"kaonTofNSigmaMax", 3.f, "Maximum absolute kaon TOF n-sigma"};
    Configurable<float> kaonCombinedNSigmaMax{"kaonCombinedNSigmaMax", 3.f, "Maximum circular TPC-TOF kaon n-sigma"};
  } pidCuts;

  struct : ConfigurableGroup {
    Configurable<int> numberOfMixedEvents{"numberOfMixedEvents", 5, "Number of mixed-event neighbours"};
    ConfigurableAxis vertexZBins{"vertexZBins", {VARIABLE_WIDTH, -10.f, -5.f, -2.f, 0.f, 2.f, 5.f, 10.f}, "Mixing bins in collision z vertex"};
    ConfigurableAxis multiplicityBins{"multiplicityBins", {VARIABLE_WIDTH, 0.f, 10.f, 20.f, 40.f, 80.f, 150.f, 300.f, 600.f, 2000.f}, "Mixing bins in FT0C multiplicity"};
    Configurable<float> nominalJetPtDifferenceMax{"nominalJetPtDifferenceMax", 2.5f, "Maximum jet-pT difference for nominal mixing (GeV/c)"};
    Configurable<int> nominalJetNchDifferenceMax{"nominalJetNchDifferenceMax", 3, "Maximum jet-constituent multiplicity difference for nominal mixing"};
    Configurable<float> nominalJetEtaDifferenceMax{"nominalJetEtaDifferenceMax", 0.1f, "Maximum jet-eta difference for nominal mixing"};
    Configurable<float> looseJetPtDifferenceMax{"looseJetPtDifferenceMax", 7.5f, "Maximum jet-pT difference for loose mixing (GeV/c)"};
    Configurable<int> looseJetNchDifferenceMax{"looseJetNchDifferenceMax", 8, "Maximum jet-constituent multiplicity difference for loose mixing"};
    Configurable<float> looseJetEtaDifferenceMax{"looseJetEtaDifferenceMax", 0.3f, "Maximum jet-eta difference for loose mixing"};
  } mixing;

  struct : ConfigurableGroup {
    Configurable<bool> fillQa{"fillQa", true, "Fill event, jet, track, and PID QA histograms"};
    Configurable<bool> fillReferenceVariations{"fillReferenceVariations", true, "Fill loose and unmatched mixed-event references"};
    Configurable<bool> fillChargeSeparated{"fillChargeSeparated", true, "Fill charge-resolved pair histograms"};
    Configurable<bool> fillPionControl{"fillPionControl", false, "Fill same-species pion control correlations"};
    Configurable<bool> fillKaonControl{"fillKaonControl", false, "Fill same-species kaon control correlations"};
  } output;

  struct : ConfigurableGroup {
    Configurable<float> nearSideDeltaPhiMax{"nearSideDeltaPhiMax", 0.2f, "Maximum absolute delta phi for near-side pair yields"};
    Configurable<float> nearSideDeltaEtaMax{"nearSideDeltaEtaMax", 0.6f, "Maximum absolute delta eta for near-side pair yields"};
    ConfigurableAxis jetPt{"jetPt", {VARIABLE_WIDTH, 0.f, 5.f, 10.f, 15.f, 20.f, 25.f, 30.f, 40.f, 50.f, 60.f, 80.f, 100.f, 120.f, 140.f}, "#it{p}_{T, jet} (GeV/#it{c})"};
    ConfigurableAxis jetEta{"jetEta", {40, -1.f, 1.f}, "#eta_{jet}"};
    ConfigurableAxis jetPhi{"jetPhi", {72, 0.f, o2::constants::math::TwoPI}, "#varphi_{jet}"};
    ConfigurableAxis jetNch{"jetNch", {61, -0.5f, 60.5f}, "#it{N}_{ch}^{jet}"};
    ConfigurableAxis baryonPt{"baryonPt", {50, 0.f, 5.f}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis baryonEta{"baryonEta", {48, -1.2f, 1.2f}, "#eta"};
    ConfigurableAxis baryonRapidity{"baryonRapidity", {48, -1.2f, 1.2f}, "#it{y}"};
    ConfigurableAxis nSigma{"nSigma", {100, -5.f, 5.f}, "n#sigma"};
    ConfigurableAxis dca{"dca", {100, -0.5f, 0.5f}, "DCA (cm)"};
    ConfigurableAxis deltaEta{"deltaEta", {80, -0.8f, 0.8f}, "#Delta#eta"};
    ConfigurableAxis deltaPhi{"deltaPhi", {80, -0.8f, 0.8f}, "#Delta#varphi"};
    ConfigurableAxis deltaR{"deltaR", {40, 0.f, 0.8f}, "#Delta#it{R}"};
    ConfigurableAxis particleJetDistance{"particleJetDistance", {40, 0.f, 0.4f}, "#Delta#it{R}(p, particle-level jet)"};
    ConfigurableAxis inclusiveDeltaEta{"inclusiveDeltaEta", {64, -1.6f, 1.6f}, "#Delta#eta"};
    ConfigurableAxis inclusiveDeltaPhi{"inclusiveDeltaPhi", {72, -o2::constants::math::PIHalf, 3.f * o2::constants::math::PIHalf}, "#Delta#varphi"};
  } axes;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::vector<int> eventSelectionBits;
  o2::aod::rctsel::RCTFlagsChecker inclusiveRctChecker;
  SliceCache mixingCache;
  SliceCache mcdMixingCache;
  SliceCache mcpMixingCache;
  SliceCache inclusiveMixingCache;
  using MixingBinning = ColumnBinningPolicy<aod::jcollision::PosZ, aod::jcollision::MultFT0C>;
  using MCPMixingBinning = ColumnBinningPolicy<aod::jmccollision::PosZ, aod::jmccollision::MultFT0C>;
  using InclusiveSameEventCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using InclusiveMixedEventCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults>;
  using InclusiveMixingBinning = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0C>;

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventCuts.eventSelections));
    inclusiveRctChecker.init(static_cast<std::string>(eventCuts.rctLabel));

    AxisSpec axisPairType{2, -0.5, 1.5, "pair type"};
    AxisSpec axisChargePairType{3, -0.5, 2.5, "charge pair type"};
    AxisSpec axisCounter{8, 0.5, 8.5, "counter"};
    AxisSpec axisMultiplicity{11, -0.5, 10.5, "multiplicity per jet"};
    AxisSpec axisJetPt{axes.jetPt};
    AxisSpec axisJetEta{axes.jetEta};
    AxisSpec axisJetPhi{axes.jetPhi};
    AxisSpec axisJetNch{axes.jetNch};
    AxisSpec axisBaryonPt{axes.baryonPt};
    AxisSpec axisBaryonEta{axes.baryonEta};
    AxisSpec axisBaryonRapidity{axes.baryonRapidity};
    AxisSpec axisNSigma{axes.nSigma};
    AxisSpec axisDca{axes.dca};
    AxisSpec axisDeltaEta{axes.deltaEta};
    AxisSpec axisDeltaPhi{axes.deltaPhi};
    AxisSpec axisDeltaR{axes.deltaR};
    AxisSpec axisParticleJetDistance{axes.particleJetDistance};
    AxisSpec axisInclusiveDeltaEta{axes.inclusiveDeltaEta};
    AxisSpec axisInclusiveDeltaPhi{axes.inclusiveDeltaPhi};
    AxisSpec axisOrigin{3, 0.5, 3.5, "origin"};

    registry.add("event/cutFlow", "Event selection;selection step;entries", HistType::kTH1F, {axisCounter});
    if (output.fillQa.value) {
      registry.add("event/zVertex", "Selected collision z vertex;#it{z}_{vtx} (cm);entries", HistType::kTH1F, {{80, -20.f, 20.f}});
      registry.add("event/multiplicityFT0C", "Selected collision FT0C multiplicity;FT0C multiplicity;entries", HistType::kTH1F, {{200, 0.f, 2000.f}});

      registry.add("jet/ptBeforePtSelection", "Charged jets before the jet-pT selection;#it{p}_{T,jet} (GeV/#it{c});entries", HistType::kTH1F, {axisJetPt});
      registry.add("jet/pt", "Selected charged jets;#it{p}_{T,jet} (GeV/#it{c});entries", HistType::kTH1F, {axisJetPt});
      registry.add("jet/eta", "Selected charged jets;#eta_{jet};entries", HistType::kTH1F, {axisJetEta});
      registry.add("jet/phi", "Selected charged jets;#varphi_{jet};entries", HistType::kTH1F, {axisJetPhi});
      registry.add("jet/nConstituents", "Charged-jet constituent multiplicity;#it{N}_{ch}^{jet};entries", HistType::kTH1F, {axisJetNch});
      registry.add("jet/nSelectedBaryons", "Selected baryon multiplicity per jet;#it{N}_{p+#bar{p}};jets", HistType::kTH1F, {axisMultiplicity});
      registry.add("jet/nProtons", "Selected proton multiplicity per jet;#it{N}_{p};jets", HistType::kTH1F, {axisMultiplicity});
      registry.add("jet/nAntiprotons", "Selected antiproton multiplicity per jet;#it{N}_{#bar{p}};jets", HistType::kTH1F, {axisMultiplicity});
      registry.add("jet/nProtonsVsAntiprotons", "Proton and antiproton multiplicities per jet;#it{N}_{p};#it{N}_{#bar{p}}", HistType::kTH2F, {axisMultiplicity, axisMultiplicity});
      registry.add("jet/nConstituentsVsPt", "Constituent multiplicity versus jet momentum;#it{p}_{T,jet} (GeV/#it{c});#it{N}_{ch}^{jet}", HistType::kTH2F, {axisJetPt, axisJetNch});
      registry.add("jet/nBaryonsVsPt", "Selected baryon multiplicity versus jet momentum;#it{p}_{T,jet} (GeV/#it{c});#it{N}_{p+#bar{p}}", HistType::kTH2F, {axisJetPt, axisMultiplicity});
      registry.add("jet/sameSignPairsVsPt", "Same-sign pair count versus jet momentum;#it{p}_{T,jet} (GeV/#it{c});pairs", HistType::kTH1F, {axisJetPt});
      registry.add("jet/oppositeSignPairsVsPt", "Opposite-sign pair count versus jet momentum;#it{p}_{T,jet} (GeV/#it{c});pairs", HistType::kTH1F, {axisJetPt});

      registry.add("proton/pt", "Selected proton and antiproton candidates;#it{p}_{T} (GeV/#it{c});entries", HistType::kTH1F, {axisBaryonPt});
      registry.add("proton/eta", "Selected proton and antiproton candidates;#eta;entries", HistType::kTH1F, {axisBaryonEta});
      registry.add("proton/rapidity", "Selected proton and antiproton candidates;#it{y};entries", HistType::kTH1F, {axisBaryonRapidity});
      registry.add("proton/tpcNSigmaVsPt", "Proton TPC PID;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{p}", HistType::kTH2F, {axisBaryonPt, axisNSigma});
      registry.add("proton/tofNSigmaVsPt", "Proton TOF PID;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{p}", HistType::kTH2F, {axisBaryonPt, axisNSigma});
      registry.add("proton/dcaXYVsPt", "Selected proton transverse DCA;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axisBaryonPt, axisDca});
      registry.add("proton/dcaZVsPt", "Selected proton longitudinal DCA;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axisBaryonPt, axisDca});
      registry.add("pid/before/tpcSignalVsSignedP", "TPC response before proton PID;charge #times #it{p} (GeV/#it{c});TPC signal", HistType::kTH2F, {{120, -6.f, 6.f}, {200, 0.f, 200.f}});
      registry.add("pid/before/tofBetaVsSignedP", "TOF response before proton PID;charge #times #it{p} (GeV/#it{c});#beta", HistType::kTH2F, {{120, -6.f, 6.f}, {120, 0.f, 1.2f}});
      registry.add("pid/before/tpcNSigmaPiVsPt", "Pion TPC PID before proton PID;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{#pi}", HistType::kTH2F, {axisBaryonPt, axisNSigma});
      registry.add("pid/before/tpcNSigmaKaVsPt", "Kaon TPC PID before proton PID;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{K}", HistType::kTH2F, {axisBaryonPt, axisNSigma});
      registry.add("pid/before/tpcNSigmaPrVsPt", "Proton TPC PID before proton PID;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}^{p}", HistType::kTH2F, {axisBaryonPt, axisNSigma});
      registry.add("pid/before/tofNSigmaPiVsPt", "Pion TOF PID before proton PID;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{#pi}", HistType::kTH2F, {axisBaryonPt, axisNSigma});
      registry.add("pid/before/tofNSigmaKaVsPt", "Kaon TOF PID before proton PID;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{K}", HistType::kTH2F, {axisBaryonPt, axisNSigma});
      registry.add("pid/before/tofNSigmaPrVsPt", "Proton TOF PID before proton PID;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}^{p}", HistType::kTH2F, {axisBaryonPt, axisNSigma});
      registry.add("pid/tofMatchDenominator", "Tracks entering PID QA;#it{p}_{T} (GeV/#it{c});tracks", HistType::kTH1F, {axisBaryonPt});
      registry.add("pid/tofMatchNumerator", "Tracks with TOF;#it{p}_{T} (GeV/#it{c});tracks", HistType::kTH1F, {axisBaryonPt});
      registry.add("pid/selectedPtByCharge", "Selected proton candidates;charge;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{2, -1.5f, 1.5f}, axisBaryonPt});
    }

    addPairHistograms("sameEvent", axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    addPairHistograms("mixedEvent/nominal", axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    addPairHistograms("mixedEvent/loose", axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    addPairHistograms("mixedEvent/random", axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    registry.add("sameEvent/chargePairSparse", "Charge-separated baryon pairs;charge pair type;#it{p}_{T,jet};#Delta#eta;#Delta#varphi;#Delta#it{R}", HistType::kTHnSparseF, {axisChargePairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR});
    registry.add("mixedEvent/nominal/chargePairSparse", "Charge-separated nominal mixed pairs;charge pair type;#it{p}_{T,jet};#Delta#eta;#Delta#varphi;#Delta#it{R}", HistType::kTHnSparseF, {axisChargePairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR});
    registry.add("mixedEvent/loose/chargePairSparse", "Charge-separated loose mixed pairs;charge pair type;#it{p}_{T,jet};#Delta#eta;#Delta#varphi;#Delta#it{R}", HistType::kTHnSparseF, {axisChargePairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR});
    registry.add("mixedEvent/random/chargePairSparse", "Charge-separated random mixed pairs;charge pair type;#it{p}_{T,jet};#Delta#eta;#Delta#varphi;#Delta#it{R}", HistType::kTHnSparseF, {axisChargePairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR});

    addInclusivePairHistograms("inclusive/sameEvent", axisChargePairType, axisInclusiveDeltaEta, axisInclusiveDeltaPhi);
    addInclusivePairHistograms("inclusive/mixedEvent", axisChargePairType, axisInclusiveDeltaEta, axisInclusiveDeltaPhi);
    registry.add("inclusive/event/multiplicity", "Inclusive selected-event multiplicity;FT0C multiplicity;events", HistType::kTH1F, {{200, 0.f, 2000.f}});
    registry.add("inclusive/normalization/sameEvent", "Inclusive same-event normalization;counter;entries", HistType::kTH1F, {axisCounter});
    registry.add("inclusive/normalization/mixedEvent", "Inclusive mixed-event normalization;counter;weighted entries", HistType::kTH1F, {axisCounter});

    registry.add("normalization/sameEvent", "Same-event normalization counters;counter;weighted entries", HistType::kTH1F, {axisCounter});
    registry.add("normalization/mixedEventNominal", "Nominal mixed-event normalization counters;counter;weighted entries", HistType::kTH1F, {axisCounter});
    registry.add("normalization/mixedEventLoose", "Loose mixed-event normalization counters;counter;weighted entries", HistType::kTH1F, {axisCounter});
    registry.add("normalization/mixedEventRandom", "Random mixed-event normalization counters;counter;weighted entries", HistType::kTH1F, {axisCounter});
    registry.add("mixedEvent/jetPtDifference", "Mixed-jet matching QA;|#Delta#it{p}_{T,jet}| (GeV/#it{c});jet pairs", HistType::kTH1F, {{80, 0.f, 40.f}});
    registry.add("mixedEvent/jetEtaDifference", "Mixed-jet matching QA;|#Delta#eta_{jet}|;jet pairs", HistType::kTH1F, {{40, 0.f, 1.f}});
    registry.add("mixedEvent/jetNchDifference", "Mixed-jet matching QA;|#Delta#it{N}_{ch}^{jet}|;jet pairs", HistType::kTH1F, {{61, -0.5f, 60.5f}});

    registry.add("mc/jet/allDetectorPt", "Accepted detector-level charged jets;#it{p}_{T,jet}^{det} (GeV/#it{c});jets", HistType::kTH1F, {axisJetPt});
    registry.add("mc/jet/allParticlePt", "Accepted particle-level charged jets;#it{p}_{T,jet}^{part} (GeV/#it{c});jets", HistType::kTH1F, {axisJetPt});
    registry.add("mc/jet/matchedDetectorPt", "Matched detector-level charged jets;#it{p}_{T,jet}^{det} (GeV/#it{c});jets", HistType::kTH1F, {axisJetPt});
    registry.add("mc/jet/matchedParticlePt", "Matched particle-level charged jets;#it{p}_{T,jet}^{part} (GeV/#it{c});jets", HistType::kTH1F, {axisJetPt});
    registry.add("mc/jet/fakeDetectorPt", "Unmatched detector-level charged jets;#it{p}_{T,jet}^{det} (GeV/#it{c});jets", HistType::kTH1F, {axisJetPt});
    registry.add("mc/jet/responsePt", "Charged-jet response;#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,jet}^{det} (GeV/#it{c})", HistType::kTH2F, {axisJetPt, axisJetPt});
    registry.add("mc/jet/relativePtResolution", "Charged-jet momentum resolution;#it{p}_{T,jet}^{part} (GeV/#it{c});(#it{p}_{T}^{det}-#it{p}_{T}^{part})/#it{p}_{T}^{part}", HistType::kTH2F, {axisJetPt, {120, -1.5f, 1.5f}});
    registry.add("mc/jet/matchingDeltaR", "Detector-particle jet matching distance;#Delta#it{R}_{match};matched jets", HistType::kTH1F, {{100, 0.f, 1.f}});
    registry.add("mc/pid/generatedPrimaryProtonPt", "Generated primary protons in accepted particle-level jets;#it{p}_{T}^{part} (GeV/#it{c});particles", HistType::kTH1F, {axisBaryonPt});
    registry.add("mc/pid/reconstructedTrueProtonPt", "Reconstructed primary true protons before PID;#it{p}_{T}^{rec} (GeV/#it{c});tracks", HistType::kTH1F, {axisBaryonPt});
    registry.add("mc/pid/selectedTrueProtonPt", "PID-selected primary true protons;#it{p}_{T}^{rec} (GeV/#it{c});tracks", HistType::kTH1F, {axisBaryonPt});
    registry.add("mc/pid/selectedCandidateOriginVsPt", "Selected proton-candidate origin;#it{p}_{T}^{rec} (GeV/#it{c});origin", HistType::kTH2F, {axisBaryonPt, axisOrigin});
    registry.add("mc/pid/selectedCandidateTruePdgVsPt", "Selected proton-candidate contamination;#it{p}_{T}^{rec} (GeV/#it{c});true PDG code", HistType::kTH2F, {axisBaryonPt, {801, -4005.f, 4005.f}});
    registry.add("mc/pid/ptResponsePrimaryProton", "Primary-proton momentum response;#it{p}_{T}^{part} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", HistType::kTH2F, {axisBaryonPt, axisBaryonPt});
    AxisSpec const axisCharge{2, -1.5f, 1.5f, "charge"};
    registry.add("mc/efficiency/generatedPrimaryProton", "Generated primary protons in matched particle-level jets;charge;#it{p}_{T}^{truth} (GeV/#it{c});#it{p}_{T,jet}^{det} (GeV/#it{c});#Delta#it{R}(p,jet_{part})", HistType::kTHnSparseF, {axisCharge, axisBaryonPt, axisJetPt, axisParticleJetDistance});
    registry.add("mc/efficiency/trackSelectedPrimaryProton", "Track-selected primary protons in matched jets;charge;#it{p}_{T}^{truth} (GeV/#it{c});#it{p}_{T,jet}^{det} (GeV/#it{c});#Delta#it{R}(p,jet_{part})", HistType::kTHnSparseF, {axisCharge, axisBaryonPt, axisJetPt, axisParticleJetDistance});
    registry.add("mc/efficiency/tofMatchedPrimaryProton", "TOF-matched primary protons in matched jets;charge;#it{p}_{T}^{truth} (GeV/#it{c});#it{p}_{T,jet}^{det} (GeV/#it{c});#Delta#it{R}(p,jet_{part})", HistType::kTHnSparseF, {axisCharge, axisBaryonPt, axisJetPt, axisParticleJetDistance});
    registry.add("mc/efficiency/pidSelectedPrimaryProton", "PID-selected primary protons in matched jets;charge;#it{p}_{T}^{truth} (GeV/#it{c});#it{p}_{T,jet}^{det} (GeV/#it{c});#Delta#it{R}(p,jet_{part})", HistType::kTHnSparseF, {axisCharge, axisBaryonPt, axisJetPt, axisParticleJetDistance});
    registry.add("mc/efficiency/generatedPair", "Generated primary-proton pairs in matched jets;pair type;#it{p}_{T,jet}^{det} (GeV/#it{c});#Delta#varphi^{truth};#Delta#it{R}^{truth}", HistType::kTHnSparseF, {axisPairType, axisJetPt, axisDeltaPhi, axisDeltaR});
    registry.add("mc/efficiency/trackSelectedPair", "Track-selected primary-proton pairs in matched jets;pair type;#it{p}_{T,jet}^{det} (GeV/#it{c});#Delta#varphi^{truth};#Delta#it{R}^{truth}", HistType::kTHnSparseF, {axisPairType, axisJetPt, axisDeltaPhi, axisDeltaR});
    registry.add("mc/efficiency/pidSelectedPair", "PID-selected primary-proton pairs in matched jets;pair type;#it{p}_{T,jet}^{det} (GeV/#it{c});#Delta#varphi^{truth};#Delta#it{R}^{truth}", HistType::kTHnSparseF, {axisPairType, axisJetPt, axisDeltaPhi, axisDeltaR});
    addPairHistograms("mc/particle/sameEvent", axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    addPairHistograms("mc/particle/mixedEvent", axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    registry.add("mc/particle/normalizationSameEvent", "Particle-level same-event normalization;counter;weighted entries", HistType::kTH1F, {axisCounter});
    registry.add("mc/particle/normalizationMixedEvent", "Particle-level mixed-event normalization;counter;weighted entries", HistType::kTH1F, {axisCounter});

    if (output.fillPionControl.value) {
      addControlSpeciesHistograms("control/pion", "Pion", axisPairType, axisCounter, axisJetPt, axisBaryonPt, axisNSigma, axisOrigin, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    }
    if (output.fillKaonControl.value) {
      addControlSpeciesHistograms("control/kaon", "Kaon", axisPairType, axisCounter, axisJetPt, axisBaryonPt, axisNSigma, axisOrigin, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    }

    setHistogramLabels();
  }

  void addControlSpeciesHistograms(std::string const& directory,
                                   std::string const& speciesName,
                                   AxisSpec const& axisPairType,
                                   AxisSpec const& axisCounter,
                                   AxisSpec const& axisJetPt,
                                   AxisSpec const& axisParticlePt,
                                   AxisSpec const& axisNSigma,
                                   AxisSpec const& axisOrigin,
                                   AxisSpec const& axisDeltaEta,
                                   AxisSpec const& axisDeltaPhi,
                                   AxisSpec const& axisDeltaR)
  {
    auto addHistogram = [&](std::string const& suffix, std::string const& title, HistType type, std::vector<AxisSpec> const& histogramAxes) {
      registry.add(directory + suffix, title.c_str(), type, histogramAxes);
    };
    AxisSpec const axisPdg{801, -4005.f, 4005.f, "true PDG code"};
    addPairHistograms(directory + "/sameEvent", axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    addPairHistograms(directory + "/mixedEvent", axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    addHistogram("/normalizationSameEvent", speciesName + " same-event normalization;counter;weighted entries", HistType::kTH1F, {axisCounter});
    addHistogram("/normalizationMixedEvent", speciesName + " mixed-event normalization;counter;weighted entries", HistType::kTH1F, {axisCounter});
    if (output.fillQa.value) {
      addHistogram("/qa/pt", "Selected " + speciesName + " candidates;#it{p}_{T} (GeV/#it{c});entries", HistType::kTH1F, {axisParticlePt});
      addHistogram("/qa/tpcNSigmaVsPt", speciesName + " TPC PID;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axisParticlePt, axisNSigma});
      addHistogram("/qa/tofNSigmaVsPt", speciesName + " TOF PID;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axisParticlePt, axisNSigma});
    }
    addHistogram("/mc/generatedPrimaryPt", "Generated primary " + speciesName + " candidates;#it{p}_{T}^{part} (GeV/#it{c});particles", HistType::kTH1F, {axisParticlePt});
    addHistogram("/mc/reconstructedTruePt", "Reconstructed true primary " + speciesName + " candidates;#it{p}_{T}^{rec} (GeV/#it{c});tracks", HistType::kTH1F, {axisParticlePt});
    addHistogram("/mc/selectedTruePt", "PID-selected true primary " + speciesName + " candidates;#it{p}_{T}^{rec} (GeV/#it{c});tracks", HistType::kTH1F, {axisParticlePt});
    addHistogram("/mc/selectedCandidateOriginVsPt", "Selected " + speciesName + " candidate origin;#it{p}_{T}^{rec} (GeV/#it{c});origin", HistType::kTH2F, {axisParticlePt, axisOrigin});
    addHistogram("/mc/selectedCandidateTruePdgVsPt", "Selected " + speciesName + " candidate contamination;#it{p}_{T}^{rec} (GeV/#it{c});true PDG code", HistType::kTH2F, {axisParticlePt, axisPdg});
    addHistogram("/mc/ptResponsePrimary", "Primary " + speciesName + " momentum response;#it{p}_{T}^{part} (GeV/#it{c});#it{p}_{T}^{rec} (GeV/#it{c})", HistType::kTH2F, {axisParticlePt, axisParticlePt});
    addPairHistograms(directory + "/mc/particle/sameEvent", axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    addPairHistograms(directory + "/mc/particle/mixedEvent", axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR);
    addHistogram("/mc/particle/normalizationSameEvent", speciesName + " particle-level same-event normalization;counter;weighted entries", HistType::kTH1F, {axisCounter});
    addHistogram("/mc/particle/normalizationMixedEvent", speciesName + " particle-level mixed-event normalization;counter;weighted entries", HistType::kTH1F, {axisCounter});
  }

  void addPairHistograms(std::string const& directory,
                         AxisSpec const& axisPairType,
                         AxisSpec const& axisJetPt,
                         AxisSpec const& axisDeltaEta,
                         AxisSpec const& axisDeltaPhi,
                         AxisSpec const& axisDeltaR)
  {
    registry.add(directory + "/sameSignDeltaEtaDeltaPhi", "Same-sign baryon pairs;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisDeltaEta, axisDeltaPhi});
    registry.add(directory + "/oppositeSignDeltaEtaDeltaPhi", "Opposite-sign baryon pairs;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisDeltaEta, axisDeltaPhi});
    registry.add(directory + "/sameSignDeltaR", "Same-sign baryon pairs;#Delta#it{R};weighted pairs", HistType::kTH1F, {axisDeltaR});
    registry.add(directory + "/oppositeSignDeltaR", "Opposite-sign baryon pairs;#Delta#it{R};weighted pairs", HistType::kTH1F, {axisDeltaR});
    registry.add(directory + "/pairSparse", "Baryon-pair distribution;pair type;#it{p}_{T,jet};#Delta#eta;#Delta#varphi;#Delta#it{R}", HistType::kTHnSparseF, {axisPairType, axisJetPt, axisDeltaEta, axisDeltaPhi, axisDeltaR});
  }

  void addInclusivePairHistograms(std::string const& directory,
                                  AxisSpec const& axisChargePairType,
                                  AxisSpec const& axisDeltaEta,
                                  AxisSpec const& axisDeltaPhi)
  {
    registry.add(directory + "/sameSignDeltaEtaDeltaPhi", "Inclusive same-sign baryon pairs;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisDeltaEta, axisDeltaPhi});
    registry.add(directory + "/oppositeSignDeltaEtaDeltaPhi", "Inclusive opposite-sign baryon pairs;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axisDeltaEta, axisDeltaPhi});
    registry.add(directory + "/chargePairSparse", "Inclusive charge-separated baryon pairs;charge pair type;#Delta#eta;#Delta#varphi", HistType::kTHnSparseF, {axisChargePairType, axisDeltaEta, axisDeltaPhi});
  }

  void setHistogramLabels()
  {
    auto eventCutFlow = registry.get<TH1>(HIST("event/cutFlow"));
    eventCutFlow->GetXaxis()->SetBinLabel(1, "all");
    eventCutFlow->GetXaxis()->SetBinLabel(2, "z vertex");
    eventCutFlow->GetXaxis()->SetBinLabel(3, "occupancy");
    eventCutFlow->GetXaxis()->SetBinLabel(4, "event selection");

    auto setNormalizationLabels = [](std::shared_ptr<TH1> const& histogram) {
      histogram->GetXaxis()->SetBinLabel(1, "events or event pairs");
      histogram->GetXaxis()->SetBinLabel(2, "selected jets or jet pairs");
      histogram->GetXaxis()->SetBinLabel(3, "same-sign pairs");
      histogram->GetXaxis()->SetBinLabel(4, "opposite-sign pairs");
      histogram->GetXaxis()->SetBinLabel(5, "near-side same-sign pairs");
      histogram->GetXaxis()->SetBinLabel(6, "near-side opposite-sign pairs");
    };
    setNormalizationLabels(registry.get<TH1>(HIST("normalization/sameEvent")));
    setNormalizationLabels(registry.get<TH1>(HIST("normalization/mixedEventNominal")));
    setNormalizationLabels(registry.get<TH1>(HIST("normalization/mixedEventLoose")));
    setNormalizationLabels(registry.get<TH1>(HIST("normalization/mixedEventRandom")));
    if (output.fillPionControl.value) {
      setNormalizationLabels(registry.get<TH1>(HIST("control/pion/normalizationSameEvent")));
      setNormalizationLabels(registry.get<TH1>(HIST("control/pion/normalizationMixedEvent")));
      setNormalizationLabels(registry.get<TH1>(HIST("control/pion/mc/particle/normalizationSameEvent")));
      setNormalizationLabels(registry.get<TH1>(HIST("control/pion/mc/particle/normalizationMixedEvent")));
    }
    if (output.fillKaonControl.value) {
      setNormalizationLabels(registry.get<TH1>(HIST("control/kaon/normalizationSameEvent")));
      setNormalizationLabels(registry.get<TH1>(HIST("control/kaon/normalizationMixedEvent")));
      setNormalizationLabels(registry.get<TH1>(HIST("control/kaon/mc/particle/normalizationSameEvent")));
      setNormalizationLabels(registry.get<TH1>(HIST("control/kaon/mc/particle/normalizationMixedEvent")));
    }
  }

  template <typename CollisionType>
  bool isSelectedCollision(CollisionType const& collision, bool fillCutFlow)
  {
    if (fillCutFlow) {
      registry.fill(HIST("event/cutFlow"), 1.f);
    }
    if (std::abs(collision.posZ()) > eventCuts.vertexZMax.value) {
      return false;
    }
    if (fillCutFlow) {
      registry.fill(HIST("event/cutFlow"), 2.f);
    }
    if (collision.trackOccupancyInTimeRange() < eventCuts.occupancyMin.value || collision.trackOccupancyInTimeRange() > eventCuts.occupancyMax.value) {
      return false;
    }
    if (fillCutFlow) {
      registry.fill(HIST("event/cutFlow"), 3.f);
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, true, eventCuts.applyRctSelection.value, static_cast<std::string>(eventCuts.rctLabel))) {
      return false;
    }
    if (fillCutFlow) {
      registry.fill(HIST("event/cutFlow"), 4.f);
    }
    return true;
  }

  template <typename JetType>
  bool isSelectedJetGeometry(JetType const& jet) const
  {
    if (jet.r() != static_cast<int>(std::lround(jetCuts.jetRadius.value * 100.f)) || std::abs(jet.eta()) > jetCuts.jetEtaMax.value) {
      return false;
    }
    if (jetCuts.areaFractionMin.value >= 0.f) {
      float const fullArea = o2::constants::math::PI * jetCuts.jetRadius.value * jetCuts.jetRadius.value;
      if (jet.area() / fullArea < jetCuts.areaFractionMin.value) {
        return false;
      }
    }
    return true;
  }

  template <typename JetType>
  bool isSelectedJet(JetType const& jet) const
  {
    return isSelectedJetGeometry(jet) && jet.pt() >= jetCuts.jetPtMin.value && jet.pt() < jetCuts.jetPtMax.value;
  }

  template <typename TrackType>
  bool isSelectedTrack(TrackType const& track) const
  {
    if (trackCuts.requirePvContributor.value && !track.isPVContributor()) {
      return false;
    }
    if (!track.hasITS() || !track.hasTPC() || !track.isGlobalTrack()) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < trackCuts.tpcCrossedRowsMin.value || track.tpcChi2NCl() > trackCuts.tpcChi2Max.value || track.itsChi2NCl() > trackCuts.itsChi2Max.value) {
      return false;
    }
    if (std::abs(track.eta()) > trackCuts.protonEtaMax.value || track.pt() < trackCuts.protonPtMin.value || track.pt() >= trackCuts.protonPtMax.value) {
      return false;
    }
    if (std::abs(track.dcaXY()) > trackCuts.dcaXYMax.value || std::abs(track.dcaZ()) > trackCuts.dcaZMax.value) {
      return false;
    }
    return true;
  }

  template <typename TrackType>
  bool isSelectedProton(TrackType const& track) const
  {
    auto const pidMethod = static_cast<ProtonPidMethod>(pidCuts.method.value);
    float const tpcNSigma = std::abs(track.tpcNSigmaPr());
    bool isProton = false;

    if (pidMethod == ProtonPidMethod::TpcOnly || track.pt() < pidCuts.tofPtThreshold.value) {
      isProton = tpcNSigma < pidCuts.tpcNSigmaMax.value;
    } else if (track.hasTOF()) {
      float const tofNSigma = std::abs(track.tofNSigmaPr());
      if (pidMethod == ProtonPidMethod::TpcTofCircular) {
        isProton = std::hypot(tpcNSigma, tofNSigma) < pidCuts.combinedNSigmaMax.value;
      } else if (pidMethod == ProtonPidMethod::TpcTofRectangular) {
        isProton = tpcNSigma < pidCuts.tpcNSigmaMax.value && tofNSigma < pidCuts.tofNSigmaMax.value;
      }
    }

    if (!isProton || !pidCuts.rejectPionKaon.value) {
      return isProton;
    }

    bool const isPionCompatible = std::abs(track.tpcNSigmaPi()) < pidCuts.rejectionNSigma.value &&
                                  (track.pt() < pidCuts.tofPtThreshold.value || (track.hasTOF() && std::abs(track.tofNSigmaPi()) < pidCuts.rejectionNSigma.value));
    bool const isKaonCompatible = std::abs(track.tpcNSigmaKa()) < pidCuts.rejectionNSigma.value &&
                                  (track.pt() < pidCuts.tofPtThreshold.value || (track.hasTOF() && std::abs(track.tofNSigmaKa()) < pidCuts.rejectionNSigma.value));
    return !isPionCompatible && !isKaonCompatible;
  }

  bool passesControlPid(float pt, bool hasTof, float tpcNSigma, float tofNSigma, float tpcMaximum, float tofMaximum, float combinedMaximum) const
  {
    auto const pidMethod = static_cast<ProtonPidMethod>(pidCuts.method.value);
    tpcNSigma = std::abs(tpcNSigma);
    tofNSigma = std::abs(tofNSigma);
    if (pidMethod == ProtonPidMethod::TpcOnly || pt < pidCuts.tofPtThreshold.value) {
      return tpcNSigma < tpcMaximum;
    }
    if (!hasTof) {
      return false;
    }
    if (pidMethod == ProtonPidMethod::TpcTofCircular) {
      return std::hypot(tpcNSigma, tofNSigma) < combinedMaximum;
    }
    return tpcNSigma < tpcMaximum && tofNSigma < tofMaximum;
  }

  template <typename TrackType>
  bool isSelectedPion(TrackType const& track) const
  {
    return passesControlPid(track.pt(), track.hasTOF(), track.tpcNSigmaPi(), track.tofNSigmaPi(), pidCuts.pionTpcNSigmaMax.value, pidCuts.pionTofNSigmaMax.value, pidCuts.pionCombinedNSigmaMax.value);
  }

  template <typename TrackType>
  bool isSelectedKaon(TrackType const& track) const
  {
    return passesControlPid(track.pt(), track.hasTOF(), track.tpcNSigmaKa(), track.tofNSigmaKa(), pidCuts.kaonTpcNSigmaMax.value, pidCuts.kaonTofNSigmaMax.value, pidCuts.kaonCombinedNSigmaMax.value);
  }

  template <typename TrackType, typename JetType>
  SelectedParticles collectSelectedParticles(JetType const& jet, bool fillQa, float weight = 1.f, bool collectControls = true)
  {
    SelectedParticles selected;
    for (auto const& jetTrack : jet.template tracks_as<JetTracksWithOriginal>()) {
      auto track = jetTrack.template track_as<TrackType>();
      if (!isSelectedTrack(track)) {
        continue;
      }

      if (fillQa) {
        float const signedMomentum = static_cast<float>(track.sign()) * track.p();
        registry.fill(HIST("pid/before/tpcSignalVsSignedP"), signedMomentum, track.tpcSignal(), weight);
        registry.fill(HIST("pid/before/tpcNSigmaPiVsPt"), track.pt(), track.tpcNSigmaPi(), weight);
        registry.fill(HIST("pid/before/tpcNSigmaKaVsPt"), track.pt(), track.tpcNSigmaKa(), weight);
        registry.fill(HIST("pid/before/tpcNSigmaPrVsPt"), track.pt(), track.tpcNSigmaPr(), weight);
        registry.fill(HIST("pid/tofMatchDenominator"), track.pt(), weight);
        if (track.hasTOF()) {
          registry.fill(HIST("pid/tofMatchNumerator"), track.pt(), weight);
          registry.fill(HIST("pid/before/tofBetaVsSignedP"), signedMomentum, track.beta(), weight);
          registry.fill(HIST("pid/before/tofNSigmaPiVsPt"), track.pt(), track.tofNSigmaPi(), weight);
          registry.fill(HIST("pid/before/tofNSigmaKaVsPt"), track.pt(), track.tofNSigmaKa(), weight);
          registry.fill(HIST("pid/before/tofNSigmaPrVsPt"), track.pt(), track.tofNSigmaPr(), weight);
        }
      }

      SelectedParticle const particle{track.pt(), track.eta(), track.phi(), static_cast<int8_t>(track.sign())};
      if (isSelectedProton(track)) {
        selected.baryons.push_back(particle);
        if (fillQa) {
          float const rapidity = RecoDecay::y(std::array{track.px(), track.py(), track.pz()}, o2::constants::physics::MassProton);
          registry.fill(HIST("proton/pt"), particle.pt, weight);
          registry.fill(HIST("proton/eta"), particle.eta, weight);
          registry.fill(HIST("proton/rapidity"), rapidity, weight);
          registry.fill(HIST("proton/tpcNSigmaVsPt"), track.pt(), track.tpcNSigmaPr(), weight);
          if (track.hasTOF()) {
            registry.fill(HIST("proton/tofNSigmaVsPt"), track.pt(), track.tofNSigmaPr(), weight);
          }
          registry.fill(HIST("proton/dcaXYVsPt"), track.pt(), track.dcaXY(), weight);
          registry.fill(HIST("proton/dcaZVsPt"), track.pt(), track.dcaZ(), weight);
          registry.fill(HIST("pid/selectedPtByCharge"), static_cast<float>(particle.sign), particle.pt, weight);
        }
      }
      if (collectControls && output.fillPionControl.value && isSelectedPion(track)) {
        selected.pions.push_back(particle);
        if (fillQa) {
          registry.fill(HIST("control/pion/qa/pt"), particle.pt, weight);
          registry.fill(HIST("control/pion/qa/tpcNSigmaVsPt"), particle.pt, track.tpcNSigmaPi(), weight);
          if (track.hasTOF()) {
            registry.fill(HIST("control/pion/qa/tofNSigmaVsPt"), particle.pt, track.tofNSigmaPi(), weight);
          }
        }
      }
      if (collectControls && output.fillKaonControl.value && isSelectedKaon(track)) {
        selected.kaons.push_back(particle);
        if (fillQa) {
          registry.fill(HIST("control/kaon/qa/pt"), particle.pt, weight);
          registry.fill(HIST("control/kaon/qa/tpcNSigmaVsPt"), particle.pt, track.tpcNSigmaKa(), weight);
          if (track.hasTOF()) {
            registry.fill(HIST("control/kaon/qa/tofNSigmaVsPt"), particle.pt, track.tofNSigmaKa(), weight);
          }
        }
      }
    }
    return selected;
  }

  static PairType getPairType(SelectedParticle const& baryon1, SelectedParticle const& baryon2)
  {
    return baryon1.sign == baryon2.sign ? PairType::SameSign : PairType::OppositeSign;
  }

  static ChargePairType getChargePairType(SelectedParticle const& baryon1, SelectedParticle const& baryon2)
  {
    if (baryon1.sign > 0 && baryon2.sign > 0) {
      return ChargePairType::PositivePositive;
    }
    if (baryon1.sign < 0 && baryon2.sign < 0) {
      return ChargePairType::NegativeNegative;
    }
    return ChargePairType::OppositeSign;
  }

  template <ControlSpecies Species>
  void fillControlNormalization(bool isMixed, float counter, float weight)
  {
    if constexpr (Species == ControlSpecies::Pion) {
      if (isMixed) {
        registry.fill(HIST("control/pion/normalizationMixedEvent"), counter, weight);
      } else {
        registry.fill(HIST("control/pion/normalizationSameEvent"), counter, weight);
      }
    } else {
      if (isMixed) {
        registry.fill(HIST("control/kaon/normalizationMixedEvent"), counter, weight);
      } else {
        registry.fill(HIST("control/kaon/normalizationSameEvent"), counter, weight);
      }
    }
  }

  template <ControlSpecies Species>
  void fillControlParticleNormalization(bool isMixed, float counter, float weight)
  {
    if constexpr (Species == ControlSpecies::Pion) {
      if (isMixed) {
        registry.fill(HIST("control/pion/mc/particle/normalizationMixedEvent"), counter, weight);
      } else {
        registry.fill(HIST("control/pion/mc/particle/normalizationSameEvent"), counter, weight);
      }
    } else {
      if (isMixed) {
        registry.fill(HIST("control/kaon/mc/particle/normalizationMixedEvent"), counter, weight);
      } else {
        registry.fill(HIST("control/kaon/mc/particle/normalizationSameEvent"), counter, weight);
      }
    }
  }

  template <ControlSpecies Species, bool IsParticleLevel>
  void fillControlPair(bool isMixed, PairType pairType, float jetPt, float deltaEta, float deltaPhi, float deltaR, float weight)
  {
    auto fillPairHistograms = [&](auto const& sameSignName, auto const& oppositeSignName, auto const& deltaRSameSignName, auto const& deltaROppositeSignName, auto const& sparseName) {
      if (pairType == PairType::SameSign) {
        registry.fill(sameSignName, deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(sameSignName, -deltaEta, -deltaPhi, 0.5f * weight);
        registry.fill(deltaRSameSignName, deltaR, weight);
      } else {
        registry.fill(oppositeSignName, deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(oppositeSignName, -deltaEta, -deltaPhi, 0.5f * weight);
        registry.fill(deltaROppositeSignName, deltaR, weight);
      }
      registry.fill(sparseName, static_cast<float>(pairType), jetPt, deltaEta, deltaPhi, deltaR, 0.5f * weight);
      registry.fill(sparseName, static_cast<float>(pairType), jetPt, -deltaEta, -deltaPhi, deltaR, 0.5f * weight);
    };

    if constexpr (Species == ControlSpecies::Pion && !IsParticleLevel) {
      if (isMixed) {
        fillPairHistograms(HIST("control/pion/mixedEvent/sameSignDeltaEtaDeltaPhi"), HIST("control/pion/mixedEvent/oppositeSignDeltaEtaDeltaPhi"), HIST("control/pion/mixedEvent/sameSignDeltaR"), HIST("control/pion/mixedEvent/oppositeSignDeltaR"), HIST("control/pion/mixedEvent/pairSparse"));
      } else {
        fillPairHistograms(HIST("control/pion/sameEvent/sameSignDeltaEtaDeltaPhi"), HIST("control/pion/sameEvent/oppositeSignDeltaEtaDeltaPhi"), HIST("control/pion/sameEvent/sameSignDeltaR"), HIST("control/pion/sameEvent/oppositeSignDeltaR"), HIST("control/pion/sameEvent/pairSparse"));
      }
    } else if constexpr (Species == ControlSpecies::Kaon && !IsParticleLevel) {
      if (isMixed) {
        fillPairHistograms(HIST("control/kaon/mixedEvent/sameSignDeltaEtaDeltaPhi"), HIST("control/kaon/mixedEvent/oppositeSignDeltaEtaDeltaPhi"), HIST("control/kaon/mixedEvent/sameSignDeltaR"), HIST("control/kaon/mixedEvent/oppositeSignDeltaR"), HIST("control/kaon/mixedEvent/pairSparse"));
      } else {
        fillPairHistograms(HIST("control/kaon/sameEvent/sameSignDeltaEtaDeltaPhi"), HIST("control/kaon/sameEvent/oppositeSignDeltaEtaDeltaPhi"), HIST("control/kaon/sameEvent/sameSignDeltaR"), HIST("control/kaon/sameEvent/oppositeSignDeltaR"), HIST("control/kaon/sameEvent/pairSparse"));
      }
    } else if constexpr (Species == ControlSpecies::Pion) {
      if (isMixed) {
        fillPairHistograms(HIST("control/pion/mc/particle/mixedEvent/sameSignDeltaEtaDeltaPhi"), HIST("control/pion/mc/particle/mixedEvent/oppositeSignDeltaEtaDeltaPhi"), HIST("control/pion/mc/particle/mixedEvent/sameSignDeltaR"), HIST("control/pion/mc/particle/mixedEvent/oppositeSignDeltaR"), HIST("control/pion/mc/particle/mixedEvent/pairSparse"));
      } else {
        fillPairHistograms(HIST("control/pion/mc/particle/sameEvent/sameSignDeltaEtaDeltaPhi"), HIST("control/pion/mc/particle/sameEvent/oppositeSignDeltaEtaDeltaPhi"), HIST("control/pion/mc/particle/sameEvent/sameSignDeltaR"), HIST("control/pion/mc/particle/sameEvent/oppositeSignDeltaR"), HIST("control/pion/mc/particle/sameEvent/pairSparse"));
      }
    } else {
      if (isMixed) {
        fillPairHistograms(HIST("control/kaon/mc/particle/mixedEvent/sameSignDeltaEtaDeltaPhi"), HIST("control/kaon/mc/particle/mixedEvent/oppositeSignDeltaEtaDeltaPhi"), HIST("control/kaon/mc/particle/mixedEvent/sameSignDeltaR"), HIST("control/kaon/mc/particle/mixedEvent/oppositeSignDeltaR"), HIST("control/kaon/mc/particle/mixedEvent/pairSparse"));
      } else {
        fillPairHistograms(HIST("control/kaon/mc/particle/sameEvent/sameSignDeltaEtaDeltaPhi"), HIST("control/kaon/mc/particle/sameEvent/oppositeSignDeltaEtaDeltaPhi"), HIST("control/kaon/mc/particle/sameEvent/sameSignDeltaR"), HIST("control/kaon/mc/particle/sameEvent/oppositeSignDeltaR"), HIST("control/kaon/mc/particle/sameEvent/pairSparse"));
      }
    }
  }

  template <ControlSpecies Species, bool IsParticleLevel>
  void fillControlSameJetPairs(std::vector<SelectedParticle> const& particles, float jetPt, float weight)
  {
    if constexpr (IsParticleLevel) {
      fillControlParticleNormalization<Species>(false, 2.f, weight);
    } else {
      fillControlNormalization<Species>(false, 2.f, weight);
    }
    for (size_t first = 0; first < particles.size(); ++first) {
      for (size_t second = first + 1; second < particles.size(); ++second) {
        PairType const pairType = getPairType(particles[first], particles[second]);
        float const deltaEta = particles[first].eta - particles[second].eta;
        float const deltaPhi = RecoDecay::constrainAngle(particles[first].phi - particles[second].phi, -o2::constants::math::PI);
        float const deltaR = std::hypot(deltaEta, deltaPhi);
        fillControlPair<Species, IsParticleLevel>(false, pairType, jetPt, deltaEta, deltaPhi, deltaR, weight);
        if constexpr (IsParticleLevel) {
          fillControlParticleNormalization<Species>(false, pairType == PairType::SameSign ? 3.f : 4.f, weight);
        } else {
          fillControlNormalization<Species>(false, pairType == PairType::SameSign ? 3.f : 4.f, weight);
        }
      }
    }
  }

  template <ControlSpecies Species, bool IsParticleLevel, typename JetType1, typename JetType2>
  void fillControlMixedJetPairs(std::vector<SelectedParticle> const& particles1,
                                std::vector<SelectedParticle> const& particles2,
                                JetType1 const& jet1,
                                JetType2 const& jet2,
                                float weight)
  {
    if constexpr (IsParticleLevel) {
      fillControlParticleNormalization<Species>(true, 2.f, weight);
    } else {
      fillControlNormalization<Species>(true, 2.f, weight);
    }
    float const mixedJetPt = 0.5f * (jet1.pt() + jet2.pt());
    for (auto const& particle1 : particles1) {
      float const localEta1 = particle1.eta - jet1.eta();
      float const localPhi1 = RecoDecay::constrainAngle(particle1.phi - jet1.phi(), -o2::constants::math::PI);
      for (auto const& particle2 : particles2) {
        float const localEta2 = particle2.eta - jet2.eta();
        float const localPhi2 = RecoDecay::constrainAngle(particle2.phi - jet2.phi(), -o2::constants::math::PI);
        float const deltaEta = localEta1 - localEta2;
        float const deltaPhi = RecoDecay::constrainAngle(localPhi1 - localPhi2, -o2::constants::math::PI);
        float const deltaR = std::hypot(deltaEta, deltaPhi);
        PairType const pairType = getPairType(particle1, particle2);
        fillControlPair<Species, IsParticleLevel>(true, pairType, mixedJetPt, deltaEta, deltaPhi, deltaR, weight);
        if constexpr (IsParticleLevel) {
          fillControlParticleNormalization<Species>(true, pairType == PairType::SameSign ? 3.f : 4.f, weight);
        } else {
          fillControlNormalization<Species>(true, pairType == PairType::SameSign ? 3.f : 4.f, weight);
        }
      }
    }
  }

  template <typename TracksType>
  std::vector<SelectedParticle> collectInclusiveBaryons(TracksType const& tracks)
  {
    std::vector<SelectedParticle> baryons;
    for (auto const& track : tracks) {
      if (!isSelectedTrack(track)) {
        continue;
      }
      if (!isSelectedProton(track)) {
        continue;
      }
      baryons.push_back({track.pt(), track.eta(), track.phi(), static_cast<int8_t>(track.sign())});
    }
    return baryons;
  }

  template <typename CollisionType>
  bool isSelectedInclusiveCollision(CollisionType const& collision)
  {
    if (std::abs(collision.posZ()) > eventCuts.vertexZMax.value ||
        collision.trackOccupancyInTimeRange() < eventCuts.occupancyMin.value ||
        collision.trackOccupancyInTimeRange() > eventCuts.occupancyMax.value) {
      return false;
    }

    uint16_t const collisionSelection = jetderiveddatautilities::setEventSelectionBit(collision);
    bool isOrCondition = false;
    bool passedOrCondition = false;
    for (auto const& selectionBit : eventSelectionBits) {
      if (selectionBit == -1) {
        isOrCondition = true;
        continue;
      }
      bool const passed = collisionSelection & (1U << selectionBit);
      if (!isOrCondition && !passed) {
        return false;
      }
      if (isOrCondition && passed) {
        passedOrCondition = true;
      }
    }
    if (isOrCondition && !passedOrCondition) {
      return false;
    }
    return !eventCuts.applyRctSelection.value || inclusiveRctChecker.checkTable(collision);
  }

  void fillInclusivePair(bool isMixed,
                         SelectedParticle const& baryon1,
                         SelectedParticle const& baryon2,
                         float weight)
  {
    PairType const pairType = getPairType(baryon1, baryon2);
    ChargePairType const chargePairType = getChargePairType(baryon1, baryon2);
    float const deltaEta = baryon1.eta - baryon2.eta;
    float const deltaPhi = RecoDecay::constrainAngle(baryon1.phi - baryon2.phi, -o2::constants::math::PIHalf);
    float const reflectedDeltaPhi = RecoDecay::constrainAngle(-deltaPhi, -o2::constants::math::PIHalf);
    if (isMixed) {
      if (pairType == PairType::SameSign) {
        registry.fill(HIST("inclusive/mixedEvent/sameSignDeltaEtaDeltaPhi"), deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(HIST("inclusive/mixedEvent/sameSignDeltaEtaDeltaPhi"), -deltaEta, reflectedDeltaPhi, 0.5f * weight);
      } else {
        registry.fill(HIST("inclusive/mixedEvent/oppositeSignDeltaEtaDeltaPhi"), deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(HIST("inclusive/mixedEvent/oppositeSignDeltaEtaDeltaPhi"), -deltaEta, reflectedDeltaPhi, 0.5f * weight);
      }
      if (output.fillChargeSeparated.value) {
        registry.fill(HIST("inclusive/mixedEvent/chargePairSparse"), static_cast<float>(chargePairType), deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(HIST("inclusive/mixedEvent/chargePairSparse"), static_cast<float>(chargePairType), -deltaEta, reflectedDeltaPhi, 0.5f * weight);
      }
    } else {
      if (pairType == PairType::SameSign) {
        registry.fill(HIST("inclusive/sameEvent/sameSignDeltaEtaDeltaPhi"), deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(HIST("inclusive/sameEvent/sameSignDeltaEtaDeltaPhi"), -deltaEta, reflectedDeltaPhi, 0.5f * weight);
      } else {
        registry.fill(HIST("inclusive/sameEvent/oppositeSignDeltaEtaDeltaPhi"), deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(HIST("inclusive/sameEvent/oppositeSignDeltaEtaDeltaPhi"), -deltaEta, reflectedDeltaPhi, 0.5f * weight);
      }
      if (output.fillChargeSeparated.value) {
        registry.fill(HIST("inclusive/sameEvent/chargePairSparse"), static_cast<float>(chargePairType), deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(HIST("inclusive/sameEvent/chargePairSparse"), static_cast<float>(chargePairType), -deltaEta, reflectedDeltaPhi, 0.5f * weight);
      }
    }
  }

  template <MixedReference Reference>
  bool isMatchedJetPair(float jetPtDifference, int jetNchDifference, float jetEtaDifference) const
  {
    if constexpr (Reference == MixedReference::Random) {
      return true;
    } else if constexpr (Reference == MixedReference::Nominal) {
      return jetPtDifference <= mixing.nominalJetPtDifferenceMax.value &&
             jetNchDifference <= mixing.nominalJetNchDifferenceMax.value &&
             jetEtaDifference <= mixing.nominalJetEtaDifferenceMax.value;
    } else {
      return jetPtDifference <= mixing.looseJetPtDifferenceMax.value &&
             jetNchDifference <= mixing.looseJetNchDifferenceMax.value &&
             jetEtaDifference <= mixing.looseJetEtaDifferenceMax.value;
    }
  }

  template <MixedReference Reference>
  void fillPair(PairType pairType, ChargePairType chargePairType, float jetPt, float deltaEta, float deltaPhi, float deltaR, float weight)
  {
    auto const pairTypeValue = static_cast<float>(pairType);
    auto fillSymmetrized = [&](auto const& sameSignName, auto const& oppositeSignName, auto const& sparseName) {
      if (pairType == PairType::SameSign) {
        registry.fill(sameSignName, deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(sameSignName, -deltaEta, -deltaPhi, 0.5f * weight);
      } else {
        registry.fill(oppositeSignName, deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(oppositeSignName, -deltaEta, -deltaPhi, 0.5f * weight);
      }
      registry.fill(sparseName, pairTypeValue, jetPt, deltaEta, deltaPhi, deltaR, 0.5f * weight);
      registry.fill(sparseName, pairTypeValue, jetPt, -deltaEta, -deltaPhi, deltaR, 0.5f * weight);
    };

    if constexpr (Reference == MixedReference::Nominal) {
      fillSymmetrized(HIST("mixedEvent/nominal/sameSignDeltaEtaDeltaPhi"), HIST("mixedEvent/nominal/oppositeSignDeltaEtaDeltaPhi"), HIST("mixedEvent/nominal/pairSparse"));
      if (pairType == PairType::SameSign) {
        registry.fill(HIST("mixedEvent/nominal/sameSignDeltaR"), deltaR, weight);
      } else {
        registry.fill(HIST("mixedEvent/nominal/oppositeSignDeltaR"), deltaR, weight);
      }
      if (output.fillChargeSeparated.value) {
        registry.fill(HIST("mixedEvent/nominal/chargePairSparse"), static_cast<float>(chargePairType), jetPt, deltaEta, deltaPhi, deltaR, 0.5f * weight);
        registry.fill(HIST("mixedEvent/nominal/chargePairSparse"), static_cast<float>(chargePairType), jetPt, -deltaEta, -deltaPhi, deltaR, 0.5f * weight);
      }
    } else if constexpr (Reference == MixedReference::Loose) {
      fillSymmetrized(HIST("mixedEvent/loose/sameSignDeltaEtaDeltaPhi"), HIST("mixedEvent/loose/oppositeSignDeltaEtaDeltaPhi"), HIST("mixedEvent/loose/pairSparse"));
      if (pairType == PairType::SameSign) {
        registry.fill(HIST("mixedEvent/loose/sameSignDeltaR"), deltaR, weight);
      } else {
        registry.fill(HIST("mixedEvent/loose/oppositeSignDeltaR"), deltaR, weight);
      }
      if (output.fillChargeSeparated.value) {
        registry.fill(HIST("mixedEvent/loose/chargePairSparse"), static_cast<float>(chargePairType), jetPt, deltaEta, deltaPhi, deltaR, 0.5f * weight);
        registry.fill(HIST("mixedEvent/loose/chargePairSparse"), static_cast<float>(chargePairType), jetPt, -deltaEta, -deltaPhi, deltaR, 0.5f * weight);
      }
    } else {
      fillSymmetrized(HIST("mixedEvent/random/sameSignDeltaEtaDeltaPhi"), HIST("mixedEvent/random/oppositeSignDeltaEtaDeltaPhi"), HIST("mixedEvent/random/pairSparse"));
      if (pairType == PairType::SameSign) {
        registry.fill(HIST("mixedEvent/random/sameSignDeltaR"), deltaR, weight);
      } else {
        registry.fill(HIST("mixedEvent/random/oppositeSignDeltaR"), deltaR, weight);
      }
      if (output.fillChargeSeparated.value) {
        registry.fill(HIST("mixedEvent/random/chargePairSparse"), static_cast<float>(chargePairType), jetPt, deltaEta, deltaPhi, deltaR, 0.5f * weight);
        registry.fill(HIST("mixedEvent/random/chargePairSparse"), static_cast<float>(chargePairType), jetPt, -deltaEta, -deltaPhi, deltaR, 0.5f * weight);
      }
    }
  }

  template <MixedReference Reference>
  void fillMixedNormalization(float counter, float weight)
  {
    if constexpr (Reference == MixedReference::Nominal) {
      registry.fill(HIST("normalization/mixedEventNominal"), counter, weight);
    } else if constexpr (Reference == MixedReference::Loose) {
      registry.fill(HIST("normalization/mixedEventLoose"), counter, weight);
    } else {
      registry.fill(HIST("normalization/mixedEventRandom"), counter, weight);
    }
  }

  void fillSameEventPair(PairType pairType, ChargePairType chargePairType, float jetPt, float deltaEta, float deltaPhi, float deltaR, float weight)
  {
    auto const pairTypeValue = static_cast<float>(pairType);
    if (pairType == PairType::SameSign) {
      registry.fill(HIST("sameEvent/sameSignDeltaEtaDeltaPhi"), deltaEta, deltaPhi, 0.5f * weight);
      registry.fill(HIST("sameEvent/sameSignDeltaEtaDeltaPhi"), -deltaEta, -deltaPhi, 0.5f * weight);
      registry.fill(HIST("sameEvent/sameSignDeltaR"), deltaR, weight);
    } else {
      registry.fill(HIST("sameEvent/oppositeSignDeltaEtaDeltaPhi"), deltaEta, deltaPhi, 0.5f * weight);
      registry.fill(HIST("sameEvent/oppositeSignDeltaEtaDeltaPhi"), -deltaEta, -deltaPhi, 0.5f * weight);
      registry.fill(HIST("sameEvent/oppositeSignDeltaR"), deltaR, weight);
    }
    registry.fill(HIST("sameEvent/pairSparse"), pairTypeValue, jetPt, deltaEta, deltaPhi, deltaR, 0.5f * weight);
    registry.fill(HIST("sameEvent/pairSparse"), pairTypeValue, jetPt, -deltaEta, -deltaPhi, deltaR, 0.5f * weight);
    if (output.fillChargeSeparated.value) {
      registry.fill(HIST("sameEvent/chargePairSparse"), static_cast<float>(chargePairType), jetPt, deltaEta, deltaPhi, deltaR, 0.5f * weight);
      registry.fill(HIST("sameEvent/chargePairSparse"), static_cast<float>(chargePairType), jetPt, -deltaEta, -deltaPhi, deltaR, 0.5f * weight);
    }
  }

  template <typename JetType>
  SelectedParticles collectParticleCandidates(JetType const& jet, bool fillQa, float weight)
  {
    SelectedParticles selected;
    for (auto const& particle : jet.template tracks_as<aod::JetParticles>()) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (particle.pt() < trackCuts.protonPtMin.value || particle.pt() >= trackCuts.protonPtMax.value ||
          std::abs(particle.eta()) > trackCuts.protonEtaMax.value) {
        continue;
      }
      int8_t const sign = particle.pdgCode() > 0 ? 1 : -1;
      SelectedParticle const selectedParticle{particle.pt(), particle.eta(), particle.phi(), sign};
      int const absolutePdg = std::abs(particle.pdgCode());
      if (absolutePdg == PDG_t::kProton) {
        selected.baryons.push_back(selectedParticle);
        if (fillQa) {
          registry.fill(HIST("mc/pid/generatedPrimaryProtonPt"), selectedParticle.pt, weight);
        }
      } else if (output.fillPionControl.value && absolutePdg == PDG_t::kPiPlus) {
        selected.pions.push_back(selectedParticle);
        if (fillQa) {
          registry.fill(HIST("control/pion/mc/generatedPrimaryPt"), selectedParticle.pt, weight);
        }
      } else if (output.fillKaonControl.value && absolutePdg == PDG_t::kKPlus) {
        selected.kaons.push_back(selectedParticle);
        if (fillQa) {
          registry.fill(HIST("control/kaon/mc/generatedPrimaryPt"), selectedParticle.pt, weight);
        }
      }
    }
    return selected;
  }

  void fillParticlePair(bool isMixed, PairType pairType, float jetPt, float deltaEta, float deltaPhi, float deltaR, float weight)
  {
    if (isMixed) {
      if (pairType == PairType::SameSign) {
        registry.fill(HIST("mc/particle/mixedEvent/sameSignDeltaEtaDeltaPhi"), deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(HIST("mc/particle/mixedEvent/sameSignDeltaEtaDeltaPhi"), -deltaEta, -deltaPhi, 0.5f * weight);
        registry.fill(HIST("mc/particle/mixedEvent/sameSignDeltaR"), deltaR, weight);
      } else {
        registry.fill(HIST("mc/particle/mixedEvent/oppositeSignDeltaEtaDeltaPhi"), deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(HIST("mc/particle/mixedEvent/oppositeSignDeltaEtaDeltaPhi"), -deltaEta, -deltaPhi, 0.5f * weight);
        registry.fill(HIST("mc/particle/mixedEvent/oppositeSignDeltaR"), deltaR, weight);
      }
      registry.fill(HIST("mc/particle/mixedEvent/pairSparse"), static_cast<float>(pairType), jetPt, deltaEta, deltaPhi, deltaR, 0.5f * weight);
      registry.fill(HIST("mc/particle/mixedEvent/pairSparse"), static_cast<float>(pairType), jetPt, -deltaEta, -deltaPhi, deltaR, 0.5f * weight);
    } else {
      if (pairType == PairType::SameSign) {
        registry.fill(HIST("mc/particle/sameEvent/sameSignDeltaEtaDeltaPhi"), deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(HIST("mc/particle/sameEvent/sameSignDeltaEtaDeltaPhi"), -deltaEta, -deltaPhi, 0.5f * weight);
        registry.fill(HIST("mc/particle/sameEvent/sameSignDeltaR"), deltaR, weight);
      } else {
        registry.fill(HIST("mc/particle/sameEvent/oppositeSignDeltaEtaDeltaPhi"), deltaEta, deltaPhi, 0.5f * weight);
        registry.fill(HIST("mc/particle/sameEvent/oppositeSignDeltaEtaDeltaPhi"), -deltaEta, -deltaPhi, 0.5f * weight);
        registry.fill(HIST("mc/particle/sameEvent/oppositeSignDeltaR"), deltaR, weight);
      }
      registry.fill(HIST("mc/particle/sameEvent/pairSparse"), static_cast<float>(pairType), jetPt, deltaEta, deltaPhi, deltaR, 0.5f * weight);
      registry.fill(HIST("mc/particle/sameEvent/pairSparse"), static_cast<float>(pairType), jetPt, -deltaEta, -deltaPhi, deltaR, 0.5f * weight);
    }
  }

  template <EfficiencyStage Stage>
  void fillEfficiencyPairs(std::vector<SelectedParticle> const& particles, float detectorJetPt, float weight)
  {
    for (size_t first = 0; first < particles.size(); ++first) {
      for (size_t second = first + 1; second < particles.size(); ++second) {
        auto const pairType = static_cast<float>(getPairType(particles[first], particles[second]));
        float const deltaEta = particles[first].eta - particles[second].eta;
        float const deltaPhi = RecoDecay::constrainAngle(particles[first].phi - particles[second].phi, -o2::constants::math::PI);
        float const deltaR = std::hypot(deltaEta, deltaPhi);
        if constexpr (Stage == EfficiencyStage::Generated) {
          registry.fill(HIST("mc/efficiency/generatedPair"), pairType, detectorJetPt, deltaPhi, deltaR, 0.5f * weight);
          registry.fill(HIST("mc/efficiency/generatedPair"), pairType, detectorJetPt, -deltaPhi, deltaR, 0.5f * weight);
        } else if constexpr (Stage == EfficiencyStage::TrackSelected) {
          registry.fill(HIST("mc/efficiency/trackSelectedPair"), pairType, detectorJetPt, deltaPhi, deltaR, 0.5f * weight);
          registry.fill(HIST("mc/efficiency/trackSelectedPair"), pairType, detectorJetPt, -deltaPhi, deltaR, 0.5f * weight);
        } else {
          registry.fill(HIST("mc/efficiency/pidSelectedPair"), pairType, detectorJetPt, deltaPhi, deltaR, 0.5f * weight);
          registry.fill(HIST("mc/efficiency/pidSelectedPair"), pairType, detectorJetPt, -deltaPhi, deltaR, 0.5f * weight);
        }
      }
    }
  }

  template <ControlSpecies Species, typename TrackType, typename ParticleType>
  void fillControlMCPidQa(TrackType const& track, ParticleType const& particle, float weight)
  {
    int const expectedPdg = Species == ControlSpecies::Pion ? PDG_t::kPiPlus : PDG_t::kKPlus;
    bool const isTruePrimary = std::abs(particle.pdgCode()) == expectedPdg && particle.isPhysicalPrimary();
    bool const isSelected = Species == ControlSpecies::Pion ? isSelectedPion(track) : isSelectedKaon(track);
    int origin = 2;
    if (particle.isPhysicalPrimary()) {
      origin = 1;
    } else if (particle.getProcess() == TMCProcess::kPDecay) {
      origin = 3;
    }

    if constexpr (Species == ControlSpecies::Pion) {
      if (isTruePrimary) {
        registry.fill(HIST("control/pion/mc/reconstructedTruePt"), track.pt(), weight);
      }
      if (!isSelected) {
        return;
      }
      registry.fill(HIST("control/pion/mc/selectedCandidateOriginVsPt"), track.pt(), static_cast<float>(origin), weight);
      registry.fill(HIST("control/pion/mc/selectedCandidateTruePdgVsPt"), track.pt(), static_cast<float>(particle.pdgCode()), weight);
      if (isTruePrimary) {
        registry.fill(HIST("control/pion/mc/selectedTruePt"), track.pt(), weight);
        registry.fill(HIST("control/pion/mc/ptResponsePrimary"), particle.pt(), track.pt(), weight);
      }
    } else {
      if (isTruePrimary) {
        registry.fill(HIST("control/kaon/mc/reconstructedTruePt"), track.pt(), weight);
      }
      if (!isSelected) {
        return;
      }
      registry.fill(HIST("control/kaon/mc/selectedCandidateOriginVsPt"), track.pt(), static_cast<float>(origin), weight);
      registry.fill(HIST("control/kaon/mc/selectedCandidateTruePdgVsPt"), track.pt(), static_cast<float>(particle.pdgCode()), weight);
      if (isTruePrimary) {
        registry.fill(HIST("control/kaon/mc/selectedTruePt"), track.pt(), weight);
        registry.fill(HIST("control/kaon/mc/ptResponsePrimary"), particle.pt(), track.pt(), weight);
      }
    }
  }

  template <typename JetType>
  void fillMCPidQa(JetType const& jet, float weight)
  {
    for (auto const& jetTrack : jet.template tracks_as<JetTracksWithOriginal>()) {
      auto track = jetTrack.template track_as<BaryonTracksMC>();
      if (!isSelectedTrack(track) || !track.has_mcParticle()) {
        continue;
      }
      auto particle = track.template mcParticle_as<aod::McParticles>();
      if (output.fillPionControl.value) {
        fillControlMCPidQa<ControlSpecies::Pion>(track, particle, weight);
      }
      if (output.fillKaonControl.value) {
        fillControlMCPidQa<ControlSpecies::Kaon>(track, particle, weight);
      }
      bool const isTruePrimaryProton = std::abs(particle.pdgCode()) == PDG_t::kProton && particle.isPhysicalPrimary();
      if (isTruePrimaryProton) {
        registry.fill(HIST("mc/pid/reconstructedTrueProtonPt"), track.pt(), weight);
      }
      if (!isSelectedProton(track)) {
        continue;
      }
      int origin = 2;
      if (particle.isPhysicalPrimary()) {
        origin = 1;
      } else if (particle.getProcess() == TMCProcess::kPDecay) {
        origin = 3;
      }
      registry.fill(HIST("mc/pid/selectedCandidateOriginVsPt"), track.pt(), static_cast<float>(origin), weight);
      registry.fill(HIST("mc/pid/selectedCandidateTruePdgVsPt"), track.pt(), static_cast<float>(particle.pdgCode()), weight);
      if (isTruePrimaryProton) {
        registry.fill(HIST("mc/pid/selectedTrueProtonPt"), track.pt(), weight);
        registry.fill(HIST("mc/pid/ptResponsePrimaryProton"), particle.pt(), track.pt(), weight);
      }
    }
  }

  template <typename TrackType, typename CollisionType, typename JetsType>
  void fillReconstructedSameEvent(CollisionType const& collision, JetsType const& jets, float weight = 1.f)
  {
    if (!isSelectedCollision(collision, true)) {
      return;
    }

    if (output.fillQa.value) {
      registry.fill(HIST("event/zVertex"), collision.posZ(), weight);
      registry.fill(HIST("event/multiplicityFT0C"), collision.multFT0C(), weight);
    }
    registry.fill(HIST("normalization/sameEvent"), 1.f, weight);
    if (output.fillPionControl.value) {
      fillControlNormalization<ControlSpecies::Pion>(false, 1.f, weight);
    }
    if (output.fillKaonControl.value) {
      fillControlNormalization<ControlSpecies::Kaon>(false, 1.f, weight);
    }

    for (auto const& jet : jets) {
      if (output.fillQa.value && isSelectedJetGeometry(jet)) {
        registry.fill(HIST("jet/ptBeforePtSelection"), jet.pt(), weight);
      }
      if (!isSelectedJet(jet)) {
        continue;
      }

      registry.fill(HIST("normalization/sameEvent"), 2.f, weight);
      if (output.fillQa.value) {
        registry.fill(HIST("jet/pt"), jet.pt(), weight);
        registry.fill(HIST("jet/eta"), jet.eta(), weight);
        registry.fill(HIST("jet/phi"), jet.phi(), weight);
        registry.fill(HIST("jet/nConstituents"), static_cast<float>(jet.tracksIds().size()), weight);
        registry.fill(HIST("jet/nConstituentsVsPt"), jet.pt(), static_cast<float>(jet.tracksIds().size()), weight);
      }

      auto selectedParticles = collectSelectedParticles<TrackType>(jet, output.fillQa.value, weight);
      auto const& baryons = selectedParticles.baryons;
      if (output.fillQa.value) {
        int numberOfProtons = 0;
        int numberOfAntiprotons = 0;
        for (auto const& baryon : baryons) {
          baryon.sign > 0 ? ++numberOfProtons : ++numberOfAntiprotons;
        }
        registry.fill(HIST("jet/nSelectedBaryons"), static_cast<float>(baryons.size()), weight);
        registry.fill(HIST("jet/nProtons"), numberOfProtons, weight);
        registry.fill(HIST("jet/nAntiprotons"), numberOfAntiprotons, weight);
        registry.fill(HIST("jet/nProtonsVsAntiprotons"), numberOfProtons, numberOfAntiprotons, weight);
        registry.fill(HIST("jet/nBaryonsVsPt"), jet.pt(), static_cast<float>(baryons.size()), weight);
      }

      for (size_t first = 0; first < baryons.size(); ++first) {
        for (size_t second = first + 1; second < baryons.size(); ++second) {
          PairType const pairType = getPairType(baryons[first], baryons[second]);
          ChargePairType const chargePairType = getChargePairType(baryons[first], baryons[second]);
          float const deltaEta = baryons[first].eta - baryons[second].eta;
          float const deltaPhi = RecoDecay::constrainAngle(baryons[first].phi - baryons[second].phi, -o2::constants::math::PI);
          float const deltaR = std::hypot(deltaEta, deltaPhi);
          fillSameEventPair(pairType, chargePairType, jet.pt(), deltaEta, deltaPhi, deltaR, weight);
          if (output.fillQa.value) {
            if (pairType == PairType::SameSign) {
              registry.fill(HIST("jet/sameSignPairsVsPt"), jet.pt(), weight);
            } else {
              registry.fill(HIST("jet/oppositeSignPairsVsPt"), jet.pt(), weight);
            }
          }

          float const counterBin = pairType == PairType::SameSign ? 3.f : 4.f;
          registry.fill(HIST("normalization/sameEvent"), counterBin, weight);
          if (std::abs(deltaPhi) < axes.nearSideDeltaPhiMax.value && std::abs(deltaEta) < axes.nearSideDeltaEtaMax.value) {
            registry.fill(HIST("normalization/sameEvent"), pairType == PairType::SameSign ? 5.f : 6.f, weight);
          }
        }
      }
      if (output.fillPionControl.value) {
        fillControlSameJetPairs<ControlSpecies::Pion, false>(selectedParticles.pions, jet.pt(), weight);
      }
      if (output.fillKaonControl.value) {
        fillControlSameJetPairs<ControlSpecies::Kaon, false>(selectedParticles.kaons, jet.pt(), weight);
      }
    }
  }

  void processSameEvent(aod::JetCollision const& collision,
                        ChargedJetsWithConstituents const& jets,
                        JetTracksWithOriginal const&,
                        BaryonTracks const&)
  {
    fillReconstructedSameEvent<BaryonTracks>(collision, jets);
  }
  PROCESS_SWITCH(BaryonPairInJets, processSameEvent, "Same-event identified baryon pairs in charged jets", true);

  template <MixedReference Reference, typename TrackType, typename JetType1, typename JetType2>
  void fillMixedJetPair(JetType1 const& jet1, JetType2 const& jet2, float weight)
  {
    float const jetPtDifference = std::abs(jet1.pt() - jet2.pt());
    int const jetNchDifference = std::abs(static_cast<int>(jet1.tracksIds().size()) - static_cast<int>(jet2.tracksIds().size()));
    float const jetEtaDifference = std::abs(jet1.eta() - jet2.eta());
    if (!isMatchedJetPair<Reference>(jetPtDifference, jetNchDifference, jetEtaDifference)) {
      return;
    }

    if constexpr (Reference == MixedReference::Random) {
      registry.fill(HIST("mixedEvent/jetPtDifference"), jetPtDifference, weight);
      registry.fill(HIST("mixedEvent/jetEtaDifference"), jetEtaDifference, weight);
      registry.fill(HIST("mixedEvent/jetNchDifference"), jetNchDifference, weight);
    }

    constexpr bool CollectControls = Reference == MixedReference::Nominal;
    auto selectedParticles1 = collectSelectedParticles<TrackType>(jet1, false, 1.f, CollectControls);
    auto selectedParticles2 = collectSelectedParticles<TrackType>(jet2, false, 1.f, CollectControls);
    auto const& baryons1 = selectedParticles1.baryons;
    auto const& baryons2 = selectedParticles2.baryons;
    float const mixedJetPt = 0.5f * (jet1.pt() + jet2.pt());
    fillMixedNormalization<Reference>(2.f, weight);

    for (auto const& baryon1 : baryons1) {
      float const localEta1 = baryon1.eta - jet1.eta();
      float const localPhi1 = RecoDecay::constrainAngle(baryon1.phi - jet1.phi(), -o2::constants::math::PI);
      for (auto const& baryon2 : baryons2) {
        float const localEta2 = baryon2.eta - jet2.eta();
        float const localPhi2 = RecoDecay::constrainAngle(baryon2.phi - jet2.phi(), -o2::constants::math::PI);
        float const deltaEta = localEta1 - localEta2;
        float const deltaPhi = RecoDecay::constrainAngle(localPhi1 - localPhi2, -o2::constants::math::PI);
        float const deltaR = std::hypot(deltaEta, deltaPhi);
        PairType const pairType = getPairType(baryon1, baryon2);
        ChargePairType const chargePairType = getChargePairType(baryon1, baryon2);
        fillPair<Reference>(pairType, chargePairType, mixedJetPt, deltaEta, deltaPhi, deltaR, weight);
        fillMixedNormalization<Reference>(pairType == PairType::SameSign ? 3.f : 4.f, weight);
        if (std::abs(deltaPhi) < axes.nearSideDeltaPhiMax.value && std::abs(deltaEta) < axes.nearSideDeltaEtaMax.value) {
          fillMixedNormalization<Reference>(pairType == PairType::SameSign ? 5.f : 6.f, weight);
        }
      }
    }
    if constexpr (Reference == MixedReference::Nominal) {
      if (output.fillPionControl.value) {
        fillControlMixedJetPairs<ControlSpecies::Pion, false>(selectedParticles1.pions, selectedParticles2.pions, jet1, jet2, weight);
      }
      if (output.fillKaonControl.value) {
        fillControlMixedJetPairs<ControlSpecies::Kaon, false>(selectedParticles1.kaons, selectedParticles2.kaons, jet1, jet2, weight);
      }
    }
  }

  void processMixedEvent(aod::JetCollisions const& collisions,
                         ChargedJetsWithConstituents const& jets,
                         JetTracksWithOriginal const&,
                         BaryonTracks const&)
  {
    MixingBinning mixingBinning{{mixing.vertexZBins, mixing.multiplicityBins}};
    auto jetsTuple = std::make_tuple(jets);
    SameKindPair<aod::JetCollisions, ChargedJetsWithConstituents, MixingBinning> mixedCollisionPairs{
      mixingBinning, mixing.numberOfMixedEvents.value, -1, collisions, jetsTuple, &mixingCache};

    for (auto iterator = mixedCollisionPairs.begin(); iterator != mixedCollisionPairs.end(); ++iterator) {
      auto& [collision1, jets1, collision2, jets2] = *iterator;
      if (!isSelectedCollision(collision1, false) || !isSelectedCollision(collision2, false)) {
        continue;
      }

      float const weight = 1.f / static_cast<float>(iterator.currentWindowNeighbours());
      registry.fill(HIST("normalization/mixedEventNominal"), 1.f, weight);
      if (output.fillPionControl.value) {
        fillControlNormalization<ControlSpecies::Pion>(true, 1.f, weight);
      }
      if (output.fillKaonControl.value) {
        fillControlNormalization<ControlSpecies::Kaon>(true, 1.f, weight);
      }
      if (output.fillReferenceVariations.value) {
        registry.fill(HIST("normalization/mixedEventLoose"), 1.f, weight);
        registry.fill(HIST("normalization/mixedEventRandom"), 1.f, weight);
      }

      for (auto const& jet1 : jets1) {
        if (!isSelectedJet(jet1)) {
          continue;
        }
        for (auto const& jet2 : jets2) {
          if (!isSelectedJet(jet2)) {
            continue;
          }
          fillMixedJetPair<MixedReference::Nominal, BaryonTracks>(jet1, jet2, weight);
          if (output.fillReferenceVariations.value) {
            fillMixedJetPair<MixedReference::Loose, BaryonTracks>(jet1, jet2, weight);
            fillMixedJetPair<MixedReference::Random, BaryonTracks>(jet1, jet2, weight);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(BaryonPairInJets, processMixedEvent, "Mixed-event identified baryon pairs in charged jets", true);

  void processMCD(aod::JetCollisionsMCD::iterator const& collision,
                  aod::JetMcCollisions const&,
                  ChargedMCDJetsWithConstituents const& jets,
                  JetTracksWithOriginal const&,
                  BaryonTracksMC const&,
                  aod::McParticles const&)
  {
    float const weight = collision.mcCollision().weight();
    fillReconstructedSameEvent<BaryonTracksMC>(collision, jets, weight);
    if (!isSelectedCollision(collision, false)) {
      return;
    }
    for (auto const& jet : jets) {
      if (isSelectedJet(jet)) {
        fillMCPidQa(jet, weight);
      }
    }
  }
  PROCESS_SWITCH(BaryonPairInJets, processMCD, "Detector-level MC reconstructed baryon pairs and PID truth QA", false);

  void processMixedEventMCD(aod::JetCollisionsMCD const& collisions,
                            aod::JetMcCollisions const&,
                            ChargedMCDJetsWithConstituents const& jets,
                            JetTracksWithOriginal const&,
                            BaryonTracksMC const&,
                            aod::McParticles const&)
  {
    MixingBinning mixingBinning{{mixing.vertexZBins, mixing.multiplicityBins}};
    auto jetsTuple = std::make_tuple(jets);
    SameKindPair<aod::JetCollisionsMCD, ChargedMCDJetsWithConstituents, MixingBinning> mixedCollisionPairs{
      mixingBinning, mixing.numberOfMixedEvents.value, -1, collisions, jetsTuple, &mcdMixingCache};

    for (auto iterator = mixedCollisionPairs.begin(); iterator != mixedCollisionPairs.end(); ++iterator) {
      auto& [collision1, jets1, collision2, jets2] = *iterator;
      if (!isSelectedCollision(collision1, false) || !isSelectedCollision(collision2, false)) {
        continue;
      }
      float const weight = collision1.mcCollision().weight() / static_cast<float>(iterator.currentWindowNeighbours());
      registry.fill(HIST("normalization/mixedEventNominal"), 1.f, weight);
      if (output.fillPionControl.value) {
        fillControlNormalization<ControlSpecies::Pion>(true, 1.f, weight);
      }
      if (output.fillKaonControl.value) {
        fillControlNormalization<ControlSpecies::Kaon>(true, 1.f, weight);
      }
      if (output.fillReferenceVariations.value) {
        registry.fill(HIST("normalization/mixedEventLoose"), 1.f, weight);
        registry.fill(HIST("normalization/mixedEventRandom"), 1.f, weight);
      }
      for (auto const& jet1 : jets1) {
        if (!isSelectedJet(jet1)) {
          continue;
        }
        for (auto const& jet2 : jets2) {
          if (!isSelectedJet(jet2)) {
            continue;
          }
          fillMixedJetPair<MixedReference::Nominal, BaryonTracksMC>(jet1, jet2, weight);
          if (output.fillReferenceVariations.value) {
            fillMixedJetPair<MixedReference::Loose, BaryonTracksMC>(jet1, jet2, weight);
            fillMixedJetPair<MixedReference::Random, BaryonTracksMC>(jet1, jet2, weight);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(BaryonPairInJets, processMixedEventMCD, "Detector-level MC mixed-event baryon pairs", false);

  void processMCP(aod::JetMcCollision const& collision,
                  ChargedMCPJetsWithConstituents const& jets,
                  aod::JetParticles const&)
  {
    if (std::abs(collision.posZ()) > eventCuts.vertexZMax.value) {
      return;
    }
    float const weight = collision.weight();
    registry.fill(HIST("mc/particle/normalizationSameEvent"), 1.f, weight);
    if (output.fillPionControl.value) {
      fillControlParticleNormalization<ControlSpecies::Pion>(false, 1.f, weight);
    }
    if (output.fillKaonControl.value) {
      fillControlParticleNormalization<ControlSpecies::Kaon>(false, 1.f, weight);
    }
    for (auto const& jet : jets) {
      if (!isSelectedJet(jet)) {
        continue;
      }
      registry.fill(HIST("mc/jet/allParticlePt"), jet.pt(), weight);
      registry.fill(HIST("mc/particle/normalizationSameEvent"), 2.f, weight);
      auto selectedParticles = collectParticleCandidates(jet, true, weight);
      auto const& baryons = selectedParticles.baryons;
      for (size_t first = 0; first < baryons.size(); ++first) {
        for (size_t second = first + 1; second < baryons.size(); ++second) {
          PairType const pairType = getPairType(baryons[first], baryons[second]);
          float const deltaEta = baryons[first].eta - baryons[second].eta;
          float const deltaPhi = RecoDecay::constrainAngle(baryons[first].phi - baryons[second].phi, -o2::constants::math::PI);
          float const deltaR = std::hypot(deltaEta, deltaPhi);
          fillParticlePair(false, pairType, jet.pt(), deltaEta, deltaPhi, deltaR, weight);
          registry.fill(HIST("mc/particle/normalizationSameEvent"), pairType == PairType::SameSign ? 3.f : 4.f, weight);
        }
      }
      if (output.fillPionControl.value) {
        fillControlSameJetPairs<ControlSpecies::Pion, true>(selectedParticles.pions, jet.pt(), weight);
      }
      if (output.fillKaonControl.value) {
        fillControlSameJetPairs<ControlSpecies::Kaon, true>(selectedParticles.kaons, jet.pt(), weight);
      }
    }
  }
  PROCESS_SWITCH(BaryonPairInJets, processMCP, "Particle-level MC baryon pairs", false);

  void processMixedEventMCP(aod::JetMcCollisions const& collisions,
                            ChargedMCPJetsWithConstituents const& jets,
                            aod::JetParticles const&)
  {
    MCPMixingBinning mixingBinning{{mixing.vertexZBins, mixing.multiplicityBins}};
    auto jetsTuple = std::make_tuple(jets);
    SameKindPair<aod::JetMcCollisions, ChargedMCPJetsWithConstituents, MCPMixingBinning> mixedCollisionPairs{
      mixingBinning, mixing.numberOfMixedEvents.value, -1, collisions, jetsTuple, &mcpMixingCache};
    for (auto iterator = mixedCollisionPairs.begin(); iterator != mixedCollisionPairs.end(); ++iterator) {
      auto& [collision1, jets1, collision2, jets2] = *iterator;
      if (std::abs(collision1.posZ()) > eventCuts.vertexZMax.value || std::abs(collision2.posZ()) > eventCuts.vertexZMax.value) {
        continue;
      }
      float const weight = collision1.weight() / static_cast<float>(iterator.currentWindowNeighbours());
      registry.fill(HIST("mc/particle/normalizationMixedEvent"), 1.f, weight);
      if (output.fillPionControl.value) {
        fillControlParticleNormalization<ControlSpecies::Pion>(true, 1.f, weight);
      }
      if (output.fillKaonControl.value) {
        fillControlParticleNormalization<ControlSpecies::Kaon>(true, 1.f, weight);
      }
      for (auto const& jet1 : jets1) {
        if (!isSelectedJet(jet1)) {
          continue;
        }
        for (auto const& jet2 : jets2) {
          if (!isSelectedJet(jet2)) {
            continue;
          }
          auto selectedParticles1 = collectParticleCandidates(jet1, false, weight);
          auto selectedParticles2 = collectParticleCandidates(jet2, false, weight);
          auto const& baryons1 = selectedParticles1.baryons;
          auto const& baryons2 = selectedParticles2.baryons;
          registry.fill(HIST("mc/particle/normalizationMixedEvent"), 2.f, weight);
          for (auto const& baryon1 : baryons1) {
            float const localEta1 = baryon1.eta - jet1.eta();
            float const localPhi1 = RecoDecay::constrainAngle(baryon1.phi - jet1.phi(), -o2::constants::math::PI);
            for (auto const& baryon2 : baryons2) {
              float const localEta2 = baryon2.eta - jet2.eta();
              float const localPhi2 = RecoDecay::constrainAngle(baryon2.phi - jet2.phi(), -o2::constants::math::PI);
              float const deltaEta = localEta1 - localEta2;
              float const deltaPhi = RecoDecay::constrainAngle(localPhi1 - localPhi2, -o2::constants::math::PI);
              float const deltaR = std::hypot(deltaEta, deltaPhi);
              PairType const pairType = getPairType(baryon1, baryon2);
              fillParticlePair(true, pairType, 0.5f * (jet1.pt() + jet2.pt()), deltaEta, deltaPhi, deltaR, weight);
              registry.fill(HIST("mc/particle/normalizationMixedEvent"), pairType == PairType::SameSign ? 3.f : 4.f, weight);
            }
          }
          if (output.fillPionControl.value) {
            fillControlMixedJetPairs<ControlSpecies::Pion, true>(selectedParticles1.pions, selectedParticles2.pions, jet1, jet2, weight);
          }
          if (output.fillKaonControl.value) {
            fillControlMixedJetPairs<ControlSpecies::Kaon, true>(selectedParticles1.kaons, selectedParticles2.kaons, jet1, jet2, weight);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(BaryonPairInJets, processMixedEventMCP, "Particle-level MC mixed-event baryon pairs", false);

  void processMatched(aod::JetCollisionsMCD::iterator const& collision,
                      aod::JetMcCollisions const&,
                      MatchedChargedMCDJets const& detectorJets,
                      MatchedChargedMCPJets const&,
                      JetTracksWithOriginalAndMC const&,
                      BaryonTracks const&,
                      aod::JetParticles const&)
  {
    if (!isSelectedCollision(collision, false) || !collision.has_mcCollision()) {
      return;
    }
    float const weight = collision.mcCollision().weight();
    for (auto const& detectorJet : detectorJets) {
      if (!isSelectedJet(detectorJet)) {
        continue;
      }
      registry.fill(HIST("mc/jet/allDetectorPt"), detectorJet.pt(), weight);
      bool hasAcceptedMatch = false;
      for (auto const& particleJet : detectorJet.template matchedJetGeo_as<MatchedChargedMCPJets>()) {
        if (!isSelectedJetGeometry(particleJet) || particleJet.mcCollisionId() != collision.mcCollisionId()) {
          continue;
        }
        hasAcceptedMatch = true;
        registry.fill(HIST("mc/jet/matchedDetectorPt"), detectorJet.pt(), weight);
        registry.fill(HIST("mc/jet/matchedParticlePt"), particleJet.pt(), weight);
        registry.fill(HIST("mc/jet/responsePt"), particleJet.pt(), detectorJet.pt(), weight);
        if (particleJet.pt() > 0.f) {
          registry.fill(HIST("mc/jet/relativePtResolution"), particleJet.pt(), (detectorJet.pt() - particleJet.pt()) / particleJet.pt(), weight);
        }
        float const deltaEta = detectorJet.eta() - particleJet.eta();
        float const deltaPhi = RecoDecay::constrainAngle(detectorJet.phi() - particleJet.phi(), -o2::constants::math::PI);
        registry.fill(HIST("mc/jet/matchingDeltaR"), std::hypot(deltaEta, deltaPhi), weight);

        std::unordered_set<int64_t> generatedParticleIds;
        std::vector<SelectedParticle> generatedProtons;
        for (auto const& particle : particleJet.template tracks_as<aod::JetParticles>()) {
          if (!particle.isPhysicalPrimary() || std::abs(particle.pdgCode()) != PDG_t::kProton) {
            continue;
          }
          if (particle.pt() < trackCuts.protonPtMin.value || particle.pt() >= trackCuts.protonPtMax.value || std::abs(particle.eta()) > trackCuts.protonEtaMax.value) {
            continue;
          }
          float const particleDeltaPhi = RecoDecay::constrainAngle(particle.phi() - particleJet.phi(), -o2::constants::math::PI);
          float const particleJetDistance = std::hypot(particle.eta() - particleJet.eta(), particleDeltaPhi);
          int8_t const sign = particle.pdgCode() > 0 ? 1 : -1;
          generatedParticleIds.insert(particle.globalIndex());
          generatedProtons.push_back({particle.pt(), particle.eta(), particle.phi(), sign});
          registry.fill(HIST("mc/efficiency/generatedPrimaryProton"), static_cast<float>(sign), particle.pt(), detectorJet.pt(), particleJetDistance, weight);
        }
        fillEfficiencyPairs<EfficiencyStage::Generated>(generatedProtons, detectorJet.pt(), weight);

        std::unordered_set<int64_t> trackSelectedParticleIds;
        std::unordered_set<int64_t> tofMatchedParticleIds;
        std::unordered_set<int64_t> pidSelectedParticleIds;
        std::vector<SelectedParticle> trackSelectedProtons;
        std::vector<SelectedParticle> pidSelectedProtons;
        for (auto const& jetTrack : detectorJet.template tracks_as<JetTracksWithOriginalAndMC>()) {
          if (!jetTrack.has_mcParticle()) {
            continue;
          }
          auto particle = jetTrack.template mcParticle_as<aod::JetParticles>();
          auto const particleId = particle.globalIndex();
          if (!generatedParticleIds.contains(particleId)) {
            continue;
          }
          auto track = jetTrack.template track_as<BaryonTracks>();
          if (!isSelectedTrack(track)) {
            continue;
          }
          float const particleDeltaPhi = RecoDecay::constrainAngle(particle.phi() - particleJet.phi(), -o2::constants::math::PI);
          float const particleJetDistance = std::hypot(particle.eta() - particleJet.eta(), particleDeltaPhi);
          int8_t const sign = particle.pdgCode() > 0 ? 1 : -1;
          SelectedParticle const selectedParticle{
            .pt = particle.pt(),
            .eta = particle.eta(),
            .phi = particle.phi(),
            .sign = sign};
          if (trackSelectedParticleIds.insert(particleId).second) {
            trackSelectedProtons.push_back(selectedParticle);
            registry.fill(HIST("mc/efficiency/trackSelectedPrimaryProton"), static_cast<float>(sign), particle.pt(), detectorJet.pt(), particleJetDistance, weight);
          }
          if (track.hasTOF() && tofMatchedParticleIds.insert(particleId).second) {
            registry.fill(HIST("mc/efficiency/tofMatchedPrimaryProton"), static_cast<float>(sign), particle.pt(), detectorJet.pt(), particleJetDistance, weight);
          }
          if (isSelectedProton(track) && pidSelectedParticleIds.insert(particleId).second) {
            pidSelectedProtons.push_back(selectedParticle);
            registry.fill(HIST("mc/efficiency/pidSelectedPrimaryProton"), static_cast<float>(sign), particle.pt(), detectorJet.pt(), particleJetDistance, weight);
          }
        }
        fillEfficiencyPairs<EfficiencyStage::TrackSelected>(trackSelectedProtons, detectorJet.pt(), weight);
        fillEfficiencyPairs<EfficiencyStage::PidSelected>(pidSelectedProtons, detectorJet.pt(), weight);
        break;
      }
      if (!hasAcceptedMatch) {
        registry.fill(HIST("mc/jet/fakeDetectorPt"), detectorJet.pt(), weight);
      }
    }
  }
  PROCESS_SWITCH(BaryonPairInJets, processMatched, "Matched charged-jet response and proton efficiency", false);

  void processInclusiveSameEvent(InclusiveSameEventCollisions::iterator const& collision, BaryonTracks const& tracks)
  {
    if (!isSelectedInclusiveCollision(collision)) {
      return;
    }
    registry.fill(HIST("inclusive/normalization/sameEvent"), 1.f);
    auto baryons = collectInclusiveBaryons(tracks);
    registry.fill(HIST("inclusive/normalization/sameEvent"), 2.f, static_cast<float>(baryons.size()));
    for (size_t first = 0; first < baryons.size(); ++first) {
      for (size_t second = first + 1; second < baryons.size(); ++second) {
        fillInclusivePair(false, baryons[first], baryons[second], 1.f);
        registry.fill(HIST("inclusive/normalization/sameEvent"), getPairType(baryons[first], baryons[second]) == PairType::SameSign ? 3.f : 4.f);
      }
    }
  }
  PROCESS_SWITCH(BaryonPairInJets, processInclusiveSameEvent, "Inclusive same-event identified baryon pairs", true);

  void processInclusiveMixedEvent(InclusiveMixedEventCollisions const& collisions, BaryonTracks const& tracks)
  {
    for (auto const& collision : collisions) {
      if (isSelectedInclusiveCollision(collision)) {
        registry.fill(HIST("inclusive/event/multiplicity"), collision.multFT0C());
      }
    }
    InclusiveMixingBinning mixingBinning{{mixing.vertexZBins, mixing.multiplicityBins}};
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<InclusiveMixedEventCollisions, BaryonTracks, InclusiveMixingBinning> mixedCollisionPairs{
      mixingBinning, mixing.numberOfMixedEvents.value, -1, collisions, tracksTuple, &inclusiveMixingCache};

    for (auto iterator = mixedCollisionPairs.begin(); iterator != mixedCollisionPairs.end(); ++iterator) {
      auto& [collision1, tracks1, collision2, tracks2] = *iterator;
      if (!isSelectedInclusiveCollision(collision1) || !isSelectedInclusiveCollision(collision2)) {
        continue;
      }
      float const weight = 1.f / static_cast<float>(iterator.currentWindowNeighbours());
      registry.fill(HIST("inclusive/normalization/mixedEvent"), 1.f, weight);
      auto baryons1 = collectInclusiveBaryons(tracks1);
      auto baryons2 = collectInclusiveBaryons(tracks2);
      registry.fill(HIST("inclusive/normalization/mixedEvent"), 2.f, weight * static_cast<float>(baryons1.size() + baryons2.size()));
      for (auto const& baryon1 : baryons1) {
        for (auto const& baryon2 : baryons2) {
          fillInclusivePair(true, baryon1, baryon2, weight);
          registry.fill(HIST("inclusive/normalization/mixedEvent"), getPairType(baryon1, baryon2) == PairType::SameSign ? 3.f : 4.f, weight);
        }
      }
    }
  }
  PROCESS_SWITCH(BaryonPairInJets, processInclusiveMixedEvent, "Inclusive mixed-event identified baryon pairs", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& configContext)
{
  return WorkflowSpec{adaptAnalysisTask<BaryonPairInJets>(configContext)};
}
