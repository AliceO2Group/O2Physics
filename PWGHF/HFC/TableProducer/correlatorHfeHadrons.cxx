// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file correlatorHfeHadrons.cxx
/// \brief Heavy Flavour electron-Hadron correaltor task - data-like, MC-reco and MC-Kine analyses.
/// \author Rashi Gupta <rashi.gupta@cern.ch>, IIT Indore
/// \author Ravindra Singh <ravindra.singh@cern.ch>, IIT Indore

#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFL/DataModel/ElectronSelectionTable.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TRandom3.h>

#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::hf_sel_electron;

const std::vector<double> zBins{VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0};
const std::vector<double> multBins{VARIABLE_WIDTH, 0., 200., 500.0, 5000.};
const std::vector<double> multBinsMcGen{VARIABLE_WIDTH, 0., 20., 50.0, 500.}; // In MCGen multiplicity is defined by counting primaries
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
const BinningType corrBinning{{zBins, multBins}, true};

using BinningTypeMcGen = ColumnBinningPolicy<aod::mccollision::PosZ, o2::aod::mult::MultMCFT0A>;

struct HfCorrelatorHfeHadrons {
  Produces<aod::HfEHadronPair> entryElectronHadronPair;
  Produces<aod::HfEHadronMcPair> entryElectronHadronPairmcGen;
  Produces<aod::HfElectron> entryElectron;
  Produces<aod::Hadron> entryHadron;
  // Configurables
  // Event Selection
  Configurable<float> zPvPosMax{"zPvPosMax", 10., "Maximum z of the primary vertex (cm)"};
  Configurable<bool> isRun3{"isRun3", true, "Data is from Run3 or Run2"};

  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "number of events mixed in ME process"};
  Configurable<float> invMassEEMax{"invMassEEMax", 0.14f, "max Invariant Mass for Photonic electron"};
  // Associated Hadron selection
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1f, "Transverse momentum range for associated hadron tracks"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8f, "Eta range  for associated hadron tracks"};
  Configurable<float> etaTrackMin{"etaTrackMin", -0.8f, "Eta range  for associated hadron tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 0.5f, "DCA XY cut"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1.0f, "DCA Z cut"};
  Configurable<bool> requireEmcal{"requireEmcal", true, "Require electron to be in EMCal"};

  // Sigma cut for non-EMCal electrons
  Configurable<float> tofNSigmaEl{"tofNSigmaEl", 3.0, "Sigma cut for electrons not in EMCal"};
  Configurable<float> tpcNsigmaElectronMin{"tpcNsigmaElectronMin", -0.5f, "min Electron TPCnsigma"};
  Configurable<float> tpcNsigmaElectronMax{"tpcNsigmaElectronMax", 3.0f, "max Electron TPCnsigma"};

  // Electron hadron correlation condition
  Configurable<bool> ptCondition{"ptCondition", true, "Electron pT should be greater than associate particle pT"};

  Configurable<float> eventFractionToAnalyze{"eventFractionToAnalyze", -1, "Fraction of events to analyze (use only for ME offline on very large samples)"};

  TRandom3 rnd{0};

  SliceCache cache;
  using TableCollisions = o2::soa::Filtered<o2::soa::Join<aod::Collisions, aod::Mults, aod::EvSels>>;
  using TableCollision = TableCollisions::iterator;
  using TableTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension>;
  using McGenTableCollisions = soa::Join<aod::McCollisions, aod::MultsExtraMC>;
  using McGenTableCollision = McGenTableCollisions::iterator;
  using McTableCollisions = o2::soa::Filtered<o2::soa::Join<TableCollisions, aod::McCollisionLabels>>;
  using McTableCollision = McTableCollisions::iterator;
  using McTableTracks = soa::Join<TableTracks, aod::McTrackLabels>;

  Filter collisionFilter = nabs(aod::collision::posZ) < zPvPosMax && aod::collision::numContrib > static_cast<uint16_t>(1);
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::HfSelEl> perCollision = aod::hf_sel_electron::collisionId;

  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "primary vertex z coordinate"};
  ConfigurableAxis binsDeltaEta{"binsDeltaEta", {30, -1.8, 1.8}, "#it{#Delta#eta}"};
  ConfigurableAxis binsDeltaPhi{"binsDeltaPhi", {32, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, "#it{#Delta#varphi}"};
  ConfigurableAxis binsPt{"binsPt", {50, 0.0, 50}, "#it{p_{T}}(GeV/#it{c})"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};
  ConfigurableAxis binsNSigma{"binsNSigma", {30, -15., 15.}, "#it{#sigma_{TPC}}"};
  ConfigurableAxis binsMass{"binsMass", {100, 0.0, 2.0}, "Mass (GeV/#it{c}^{2}); entries"};

  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    AxisSpec const axisPosZ = {binsPosZ, "Pos Z"};
    AxisSpec axisDeltaEta = {binsDeltaEta, "#Delta #eta = #eta_{Electron}- #eta_{Hadron}"};
    AxisSpec axisDeltaPhi = {binsDeltaPhi, "#Delta #varphi = #varphi_{Electron}- #varphi_{Hadron}"};
    AxisSpec axisPt = {binsPt, "#it{p_{T}}(GeV/#it{c})"};
    AxisSpec axisMass = {binsMass, "Mass (GeV/#it{c}^{2}); entries"};

    AxisSpec const axisPoolBin = {binsPoolBin, "PoolBin"};
    AxisSpec axisNSigma = {binsNSigma, "it{#sigma_{TPC}}"};

    registry.add("hZvertex", "z vertex", {HistType::kTH1D, {axisPosZ}});
    registry.add("hNevents", "No of events", {HistType::kTH1D, {{3, 1, 4}}});
    registry.add("hLikeMass", "Like mass", {HistType::kTH1D, {{axisMass}}});
    registry.add("hUnLikeMass", "unLike mass", {HistType::kTH1D, {{axisMass}}});
    registry.add("hLikeSignPt", "Like sign Momentum ", {HistType::kTH1D, {{axisPt}}});
    registry.add("hUnLikeSignPt", "UnLike sign Momentum", {HistType::kTH1D, {{axisPt}}});
    registry.add("hInclusiveEHCorrel", "Sparse for Delta phi and Delta eta Inclusive Electron with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hLSEHCorrel", "Sparse for Delta phi and Delta eta Like sign Electron pair  with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hULSEHCorrel", "Sparse for Delta phi and Delta eta  UnLike sign Electron pair with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});

    registry.add("hMCgenNonHfEHCorrel", "Sparse for Delta phi and Delta eta  for McGen Non Hf Electron  with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hMCgenInclusiveEHCorrl", "Sparse for Delta phi and Delta eta  for McGen Electron pair  with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hptElectron", "hptElectron", {HistType::kTH1D, {axisPt}});
    registry.add("hptHadron", "hptHadron", {HistType::kTH1D, {axisPt}});
    registry.add("hMCgenptHadron", "hMCgenptHadron", {HistType::kTH1D, {axisPt}});
    registry.add("hMCgenptHadronprimary", "hMCgenptHadronprimary", {HistType::kTH1D, {axisPt}});

    registry.add("hMixEventInclusiveEHCorrl", "Sparse for mix event Delta phi and Delta eta Inclusive Electron with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hMixEventLSEHCorrel", "Sparse for mix event Delta phi and Delta eta Like sign Electron pair with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hMixEventULSEHCorrel", "Sparse for mix event Delta phi and Delta eta Unlike sign Electron pair with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hMixEventMcGenInclusiveEHCorrl", "Sparse for mix event Delta phi and Delta eta Mc gen  Inclusive Electron with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hMixEventMcGenNonHfEHCorrl", "Sparse for mix event Delta phi and Delta eta Mc gen Non Hf Inclusive Electron with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hElectronBin", "Electron bin", {HistType::kTH1D, {axisPoolBin}});
    registry.add("hLSElectronBin", "Electron bin", {HistType::kTH1D, {axisPoolBin}});
    registry.add("hULSElectronBin", "Electron bin", {HistType::kTH1D, {axisPoolBin}});
    registry.add("hTracksBin", "Particles associated pool bin", {HistType::kTH1D, {axisPoolBin}});
  }

  // Associated Hadron Selection Cut
  template <typename T>
  bool selAssoHadron(T const& track)
  {
    if (!track.isGlobalTrackWoDCA()) {
      return false;
    }

    if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
      return false;
    }
    if (track.eta() < etaTrackMin || track.eta() > etaTrackMax) {
      return false;
    }
    if (track.pt() < ptTrackMin) {
      return false;
    }
    return true;
  }

  // Electron-hadron Correlation
  template <typename TracksType, typename ElectronType, typename CollisionType, typename BcType>
  void fillCorrelation(CollisionType const& collision, ElectronType const& electrons, TracksType const& tracks, BcType const&)
  {
    if (!(isRun3 ? collision.sel8() : (collision.sel7() && collision.alias_bit(kINT7)))) {
      return;
    }
    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFT0M()));
    auto bc = collision.template bc_as<BcType>();
    int gCollisionId = collision.globalIndex();
    int64_t timeStamp = bc.timestamp();

    bool skipEventTableFilling = false;
    if (eventFractionToAnalyze > 0) {
      if (rnd.Uniform(0, 1) > eventFractionToAnalyze) {
        skipEventTableFilling = true;
      }
    }

    registry.fill(HIST("hNevents"), 1);
    // Add hadron Table For Mix Event Electron Hadron correlation
    if (!skipEventTableFilling) {
      registry.fill(HIST("hZvertex"), collision.posZ());
      for (const auto& hTrack : tracks) {
        if (!selAssoHadron(hTrack)) {
          continue;
        }
        registry.fill(HIST("hTracksBin"), poolBin);
        registry.fill(HIST("hptHadron"), hTrack.pt());
        entryHadron(hTrack.phi(), hTrack.eta(), hTrack.pt(), poolBin, gCollisionId, timeStamp);
      }
    }

    //  Construct Deta Phi between electrons and hadrons

    double ptElectron = -999;
    double phiElectron = -999;
    double etaElectron = -999;

    for (const auto& eTrack : electrons) {
      ptElectron = eTrack.ptTrack();
      phiElectron = eTrack.phiTrack();
      etaElectron = eTrack.etaTrack();
      bool acceptElectron = false;

      double deltaPhi = -999;
      double deltaEta = -999;
      double ptHadron = -999;
      double etaHadron = -999;
      double phiHadron = -999;
      // EMCal electron
      if (eTrack.isEmcal() && requireEmcal) {
        acceptElectron = true;
      } else if (!eTrack.isEmcal() && !requireEmcal) {
        acceptElectron = true;
      }

      if (!acceptElectron) {
        continue; // skip electron if not passing criteria
      }

      registry.fill(HIST("hptElectron"), ptElectron);
      int nElectronLS = 0;
      int nElectronUS = 0;

      auto spanLS = eTrack.lSMassEE();
      auto spanUS = eTrack.uLSMassEE();

      if (!spanLS.empty()) {
        for (size_t i = 0; i < spanLS.size(); ++i) {
          float massLike = spanLS[i]; // <-- access i-th element

          registry.fill(HIST("hLSElectronBin"), poolBin);
          registry.fill(HIST("hLikeMass"), massLike);
          registry.fill(HIST("hLikeSignPt"), eTrack.ptTrack());

          if (massLike <= invMassEEMax) {
            ++nElectronLS;
          }
        }
      }
      if (!spanUS.empty()) {
        for (size_t i = 0; i < spanUS.size(); ++i) {
          float massUnLike = spanUS[i]; // <-- access i-th element

          registry.fill(HIST("hULSElectronBin"), poolBin);
          registry.fill(HIST("hUnLikeMass"), massUnLike);
          registry.fill(HIST("hUnLikeSignPt"), eTrack.ptTrack());

          if (massUnLike <= invMassEEMax) {
            ++nElectronUS;
          }
        }
      }

      if (!skipEventTableFilling) {
        registry.fill(HIST("hElectronBin"), poolBin);
        entryElectron(phiElectron, etaElectron, ptElectron, nElectronLS, nElectronUS, poolBin, gCollisionId, timeStamp);
      }

      for (const auto& hTrack : tracks) {
        // Apply Hadron cut
        if (!selAssoHadron(hTrack)) {
          continue;
        }
        ptHadron = hTrack.pt();
        phiHadron = hTrack.phi();
        etaHadron = hTrack.eta();
        if (hTrack.globalIndex() == eTrack.trackId()) {
          continue;
        }

        if (ptCondition && (ptElectron < ptHadron)) {
          continue;
        }

        deltaPhi = RecoDecay::constrainAngle(phiElectron - phiHadron, -o2::constants::math::PIHalf);
        deltaEta = etaElectron - etaHadron;
        registry.fill(HIST("hInclusiveEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);

        int nElHadLSCorr = 0;
        int nElHadUSCorr = 0;
        if (eTrack.nElPairLS() > 0) {
          for (int i = 0; i < eTrack.nElPairLS(); ++i) {
            ++nElHadLSCorr;
            registry.fill(HIST("hLSEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
          }
        }
        if (eTrack.nElPairUS() > 0) {
          for (int i = 0; i < eTrack.nElPairUS(); ++i) {

            registry.fill(HIST("hULSEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
            ++nElHadUSCorr;
          }
        }
        entryElectronHadronPair(deltaPhi, deltaEta, ptElectron, ptHadron, poolBin, nElHadLSCorr, nElHadUSCorr);

      } // end Hadron Track loop

    } // end Electron loop
  }

  // mix event electron-hadron correlation

  template <typename TracksType, typename ElectronType, typename CollisionType1, typename CollisionType2>
  void fillMixCorrelation(CollisionType1 const&, CollisionType2 const& c2, ElectronType const& tracks1, TracksType const& tracks2)
  {
    if (!(isRun3 ? c2.sel8() : (c2.sel7() && c2.alias_bit(kINT7)))) {
      return;
    }
    double ptElectronMix = -999;
    double phiElectronMix = -999;
    double etaElectronMix = -999;
    double deltaPhiMix = -999;
    double deltaEtaMix = -999;
    double ptHadronMix = -999;
    double etaHadronMix = -999;
    double phiHadronMix = -999;
    int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFT0M()));
    for (const auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
      if (!t1.isEmcal()) {
        continue;
      }

      ptHadronMix = t2.pt();
      ptElectronMix = t1.ptTrack();
      phiElectronMix = t1.phiTrack();
      phiHadronMix = t2.phi();
      etaElectronMix = t1.etaTrack();
      etaHadronMix = t2.eta();
      if (!selAssoHadron(t2)) {
        continue;
      }

      if (ptCondition && (ptElectronMix < ptHadronMix)) {
        continue;
      }

      deltaPhiMix = RecoDecay::constrainAngle(phiElectronMix - phiHadronMix, -o2::constants::math::PIHalf);
      deltaEtaMix = etaElectronMix - etaHadronMix;

      registry.fill(HIST("hMixEventInclusiveEHCorrl"), ptElectronMix, ptHadronMix, deltaPhiMix, deltaEtaMix);
      int nElHadLSCorr = 0;
      int nElHadUSCorr = 0;
      if (t1.nElPairLS() > 0) {
        for (int i = 0; i < t1.nElPairLS(); ++i) {

          registry.fill(HIST("hMixEventLSEHCorrel"), ptElectronMix, ptHadronMix, deltaPhiMix, deltaEtaMix);
          ++nElHadLSCorr;
        }
      }
      if (t1.nElPairUS() > 0) {
        for (int i = 0; i < t1.nElPairUS(); ++i) {

          registry.fill(HIST("hMixEventULSEHCorrel"), ptElectronMix, ptHadronMix, deltaPhiMix, deltaEtaMix);
          ++nElHadUSCorr;
        }
      }
      entryElectronHadronPair(deltaPhiMix, deltaEtaMix, ptElectronMix, ptHadronMix, poolBin, nElHadLSCorr, nElHadUSCorr);
    }
  }

  // =======  Process starts for Data, Same event ============

  void processData(TableCollision const& collision,
                   aod::HfCorrSelEl const& electrons,
                   TableTracks const& tracks,
                   aod::BCsWithTimestamps const& bcs)
  {
    fillCorrelation(collision, electrons, tracks, bcs);
  }

  PROCESS_SWITCH(HfCorrelatorHfeHadrons, processData, "Process for Data", true);

  // =======  Process starts for McRec, Same event ============

  void processMcRec(McTableCollision const& mcCollision,
                    aod::HfCorrSelEl const& mcElectrons,
                    McTableTracks const& mcTracks,
                    aod::BCsWithTimestamps const& bcs)
  {
    fillCorrelation(mcCollision, mcElectrons, mcTracks, bcs);
  }

  PROCESS_SWITCH(HfCorrelatorHfeHadrons, processMcRec, "Process MC Reco mode", false);

  void processMcGen(McGenTableCollision const& mcCollision, aod::McParticles const& mcParticles, aod::HfMcGenSelEl const& electrons)
  {

    BinningTypeMcGen const corrBinningMcGen{{zBins, multBinsMcGen}, true};
    int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), mcCollision.multMCFT0A()));

    for (const auto& particleMc : mcParticles) {
      if (particleMc.eta() < etaTrackMin || particleMc.eta() > etaTrackMax) {
        continue;
      }
      if (particleMc.pt() < ptTrackMin) {
        continue;
      }

      registry.fill(HIST("hMCgenptHadron"), particleMc.pt());
      if (particleMc.isPhysicalPrimary()) {

        registry.fill(HIST("hMCgenptHadronprimary"), particleMc.pt());
      }
    }

    double ptElectron = 0;
    double phiElectron = 0;
    double etaElectron = 0;

    for (const auto& electronMc : electrons) {
      double ptHadron = 0;
      double phiHadron = 0;
      double etaHadron = 0;
      double deltaPhi = 0;
      double deltaEta = 0;
      ptElectron = electronMc.ptTrackMc();
      phiElectron = electronMc.phiTrackMc();
      etaElectron = electronMc.etaTrackMc();
      for (const auto& particleMc : mcParticles) {
        if (particleMc.globalIndex() == electronMc.trackId()) {

          continue;
        }

        // Associated hadron Selection //////
        if (!particleMc.isPhysicalPrimary()) {
          continue;
        }

        if (particleMc.eta() < etaTrackMin || particleMc.eta() > etaTrackMax) {
          continue;
        }
        if (particleMc.pt() < ptTrackMin) {
          continue;
        }
        ptHadron = particleMc.pt();
        phiHadron = particleMc.phi();
        etaHadron = particleMc.eta();
        if (ptCondition && (ptElectron < ptHadron)) {
          return; // Apply pT condition
        }
        deltaPhi = RecoDecay::constrainAngle(phiElectron - phiHadron, -o2::constants::math::PIHalf);
        deltaEta = etaElectron - etaHadron;
        bool isNonHfeCorr = false;
        if (electronMc.isNonHfeMc()) {

          registry.fill(HIST("hMCgenNonHfEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
          isNonHfeCorr = true;
        } else {

          registry.fill(HIST("hMCgenInclusiveEHCorrl"), ptElectron, ptHadron, deltaPhi, deltaEta);
        }
        entryElectronHadronPairmcGen(deltaPhi, deltaEta, ptElectron, ptHadron, poolBin, isNonHfeCorr);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorHfeHadrons, processMcGen, "Process MC Gen  mode", false);
  // ====================== Implement Event mixing on Data ===============================

  // ====================== Implement Event mixing on Data ===================================

  void processDataMixedEvent(TableCollisions const& collision, aod::HfCorrSelEl const& electron, TableTracks const& tracks)
  {
    auto tracksTuple = std::make_tuple(electron, tracks);
    Pair<TableCollisions, aod::HfCorrSelEl, TableTracks, BinningType> const pair{corrBinning, numberEventsMixed, -1, collision, tracksTuple, &cache};

    // loop over the rows of the new table
    for (const auto& [c1, tracks1, c2, tracks2] : pair) {

      fillMixCorrelation(c1, c2, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(HfCorrelatorHfeHadrons, processDataMixedEvent, "Process Mixed Event Data", false);

  // ====================== Implement Event mixing on McRec ===================================

  void processMcRecMixedEvent(McTableCollisions const& mccollision, aod::HfCorrSelEl const& electron, McTableTracks const& mcTracks)
  {
    auto tracksTuple = std::make_tuple(electron, mcTracks);
    Pair<McTableCollisions, aod::HfCorrSelEl, McTableTracks, BinningType> const pairMcRec{corrBinning, numberEventsMixed, -1, mccollision, tracksTuple, &cache};

    // loop over the rows of the new table
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcRec) {

      fillMixCorrelation(c1, c2, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(HfCorrelatorHfeHadrons, processMcRecMixedEvent, "Process Mixed Event MC Reco mode", false);

  void processMcGenMixedEvent(McGenTableCollisions const& mcCollision, aod::HfMcGenSelEl const& electrons, aod::McParticles const& mcParticles)
  {

    BinningTypeMcGen const corrBinningMcGen{{zBins, multBinsMcGen}, true};

    auto tracksTuple = std::make_tuple(electrons, mcParticles);
    Pair<McGenTableCollisions, aod::HfMcGenSelEl, aod::McParticles, BinningTypeMcGen> const pairMcGen{corrBinningMcGen, 5, -1, mcCollision, tracksTuple, &cache};

    // loop over the rows of the new table
    double ptElectronMix = -999;
    double phiElectronMix = -999;
    double etaElectronMix = -999;
    double deltaPhiMix = -999;
    double deltaEtaMix = -999;
    double ptHadronMix = -999;
    double etaHadronMix = -999;
    double phiHadronMix = -999;
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcGen) {
      int poolBin = corrBinningMcGen.getBin(std::make_tuple(c1.posZ(), c1.multMCFT0A()));
      for (const auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        ptHadronMix = t2.pt();
        ptElectronMix = t1.ptTrackMc();
        phiElectronMix = t1.phiTrackMc();
        phiHadronMix = t2.phi();
        etaElectronMix = t1.etaTrackMc();
        etaHadronMix = t2.eta();
        if (t2.eta() < etaTrackMin || t2.eta() > etaTrackMax) {
          continue;
        }
        if (t2.pt() < ptTrackMin) {
          continue;
        }
        if (ptCondition && (ptElectronMix < ptHadronMix)) {
          continue;
        }

        deltaPhiMix = RecoDecay::constrainAngle(phiElectronMix - phiHadronMix, -o2::constants::math::PIHalf);
        deltaEtaMix = etaElectronMix - etaHadronMix;
        bool isNonHfeCorr = false;
        if (t1.isNonHfeMc()) {
          isNonHfeCorr = true;
          registry.fill(HIST("hMixEventMcGenNonHfEHCorrl"), ptElectronMix, ptHadronMix, deltaPhiMix, deltaEtaMix);
        } else {

          registry.fill(HIST("hMixEventMcGenInclusiveEHCorrl"), ptElectronMix, ptHadronMix, deltaPhiMix, deltaEtaMix);
        }

        entryElectronHadronPairmcGen(deltaPhiMix, deltaEtaMix, ptElectronMix, ptHadronMix, poolBin, isNonHfeCorr);
      }
    }
  }
  PROCESS_SWITCH(HfCorrelatorHfeHadrons, processMcGenMixedEvent, "Process Mixed Event MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCorrelatorHfeHadrons>(cfgc)};
}
