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

#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFL/DataModel/ElectronSelectionTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::hf_sel_electron;
std::vector<double> zBins{VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0};
std::vector<double> multBins{VARIABLE_WIDTH, 0., 200., 500.0, 5000.};
std::vector<double> multBinsMcGen{VARIABLE_WIDTH, 0., 20., 50.0, 500.}; // In MCGen multiplicity is defined by counting primaries
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
BinningType corrBinning{{zBins, multBins}, true};
using BinningTypeMcGen = ColumnBinningPolicy<aod::mccollision::PosZ, o2::aod::mult::MultMCFT0A>;
struct HfCorrelatorHfeHadrons {

  Produces<aod::HfEHadronPair> entryElectronHadronPair;
  Produces<aod::HfEHadronMcPair> entryElectronHadronPairmcGen;
  // Configurables
  // Event Selection
  Configurable<float> zPvPosMax{"zPvPosMax", 10., "Maximum z of the primary vertex (cm)"};
  Configurable<bool> isRun3{"isRun3", true, "Data is from Run3 or Run2"};

  // Associated Hadron selection
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1f, "Transverse momentum range for associated hadron tracks"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8f, "Eta range  for associated hadron tracks"};
  Configurable<float> etaTrackMin{"etaTrackMin", -0.8f, "Eta range  for associated hadron tracks"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 0.5f, "DCA XY cut"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1.0f, "DCA Z cut"};

  // Electron hadron correlation condition
  Configurable<bool> ptCondition{"ptCondition", true, "Electron pT should be greater than associate particle pT"};

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
  HistogramConfigSpec hCorrelSpec{HistType::kTHnSparseD, {{30, 0., 30.}, {20, 0., 20.}, {32, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {50, -1.8, 1.8}}};

  HistogramRegistry registry{
    "registry",
    {{"hInclusiveEHCorrel", "Sparse for Delta phi and Delta eta Inclusive Electron with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hLSEHCorrel", "Sparse for Delta phi and Delta eta Like sign Electron pair  with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hULSEHCorrel", "Sparse for Delta phi and Delta eta  UnLike sign Electron pair with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hMCgenNonHfEHCorrel", "Sparse for Delta phi and Delta eta Non Hf for McGen Inclusive Electron  with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hMCgenInclusiveEHCorrl", "Sparse for Delta phi and Delta eta  for McGen Electron pair  with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hptElectron", "hptElectron", {HistType::kTH1F, {{100, 0, 100}}}},

     {"hMixEventInclusiveEHCorrl", "Sparse for mix event Delta phi and Delta eta Inclusive Electron with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hMixEventLSEHCorrel", "Sparse for mix event Delta phi and Delta eta Like sign Electron pair with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hMixEventULSEHCorrel", "Sparse for mix event Delta phi and Delta eta Unlike sign Electron pair with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hMixEventMcGenInclusiveEHCorrl", "Sparse for mix event Delta phi and Delta eta Mc gen  Inclusive Electron with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hMixEventMcGenNonHfEHCorrl", "Sparse for mix event Delta phi and Delta eta Mc gen  Inclusive Electron with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec}}};

  void init(InitContext&)
  {
    registry.get<THnSparse>(HIST("hInclusiveEHCorrel"))->Sumw2();
    registry.get<THnSparse>(HIST("hLSEHCorrel"))->Sumw2();
    registry.get<THnSparse>(HIST("hULSEHCorrel"))->Sumw2();
    registry.get<THnSparse>(HIST("hMCgenInclusiveEHCorrl"))->Sumw2();
    registry.get<THnSparse>(HIST("hMCgenNonHfEHCorrel"))->Sumw2();
    registry.get<THnSparse>(HIST("hMixEventInclusiveEHCorrl"))->Sumw2();
    registry.get<THnSparse>(HIST("hMixEventLSEHCorrel"))->Sumw2();
    registry.get<THnSparse>(HIST("hMixEventULSEHCorrel"))->Sumw2();
    registry.get<THnSparse>(HIST("hMixEventMcGenInclusiveEHCorrl"))->Sumw2();
    registry.get<THnSparse>(HIST("hMixEventMcGenNonHfEHCorrl"))->Sumw2();
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
  template <typename TracksType, typename ElectronType, typename CollisionType>
  void fillCorrelation(CollisionType const& collision, ElectronType const& electron, TracksType const& tracks)
  {
    if (!(isRun3 ? collision.sel8() : (collision.sel7() && collision.alias_bit(kINT7))))
      return;
    int poolBin = corrBinning.getBin(std::make_tuple(collision.posZ(), collision.multFV0M()));

    //  Construct Deta Phi between electrons and hadrons

    double ptElectron = -999;
    double phiElectron = -999;
    double etaElectron = -999;

    for (const auto& eTrack : electron) {
      ptElectron = eTrack.ptTrack();
      phiElectron = eTrack.phiTrack();
      etaElectron = eTrack.etaTrack();

      double deltaPhi = -999;
      double deltaEta = -999;
      double ptHadron = -999;
      double etaHadron = -999;
      double phiHadron = -999;
      if (!eTrack.isEmcal()) {
        continue;
      }
      registry.fill(HIST("hptElectron"), ptElectron);

      for (const auto& hTrack : tracks) {
        if (hTrack.globalIndex() == eTrack.trackId()) {
          continue;
        }

        // Apply Hadron cut
        if (!selAssoHadron(hTrack)) {
          continue;
        }

        ptHadron = hTrack.pt();
        phiHadron = hTrack.phi();
        etaHadron = hTrack.eta();

        if (ptCondition && (ptElectron < ptHadron)) {
          continue;
        }

        deltaPhi = RecoDecay::constrainAngle(phiElectron - phiHadron, -o2::constants::math::PIHalf);
        deltaEta = etaElectron - etaHadron;
        registry.fill(HIST("hInclusiveEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
        int isLSElectroncorr = 0;
        int isULSElectroncorr = 0;
        if (eTrack.isLSElectron() > 0) {
          for (int i = 0; i < eTrack.isLSElectron(); ++i) {

            registry.fill(HIST("hLSEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
            ++isLSElectroncorr;
          }
        }
        if (eTrack.isULSElectron() > 0) {
          for (int i = 0; i < eTrack.isULSElectron(); ++i) {

            registry.fill(HIST("hULSEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
            ++isULSElectroncorr;
          }
        }
        entryElectronHadronPair(deltaPhi, deltaEta, ptElectron, ptHadron, poolBin, isLSElectroncorr, isULSElectroncorr);
      }
    }
  }

  // mix event electron-hadron correlation

  template <typename TracksType, typename ElectronType, typename CollisionType1, typename CollisionType2>
  void fillMixCorrelation(CollisionType1 const&, CollisionType2 const& c2, ElectronType const& tracks1, TracksType const& tracks2)
  {
    if (!(isRun3 ? c2.sel8() : (c2.sel7() && c2.alias_bit(kINT7))))
      return;
    double ptElectronMix = -999;
    double phiElectronMix = -999;
    double etaElectronMix = -999;
    double deltaPhiMix = -999;
    double deltaEtaMix = -999;
    double ptHadronMix = -999;
    double etaHadronMix = -999;
    double phiHadronMix = -999;
    int poolBin = corrBinning.getBin(std::make_tuple(c2.posZ(), c2.multFV0M()));
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
      int isLSElectroncorr = 0;
      int isULSElectroncorr = 0;
      if (t1.isLSElectron() > 0) {
        for (int i = 0; i < t1.isLSElectron(); ++i) {

          registry.fill(HIST("hMixEventLSEHCorrel"), ptElectronMix, ptHadronMix, deltaPhiMix, deltaEtaMix);
          ++isLSElectroncorr;
        }
      }
      if (t1.isULSElectron() > 0) {
        for (int i = 0; i < t1.isULSElectron(); ++i) {

          registry.fill(HIST("hMixEventULSEHCorrel"), ptElectronMix, ptHadronMix, deltaPhiMix, deltaEtaMix);
          ++isULSElectroncorr;
        }
      }
      entryElectronHadronPair(deltaPhiMix, deltaEtaMix, ptElectronMix, ptHadronMix, poolBin, isLSElectroncorr, isULSElectroncorr);
    }
  }

  // =======  Process starts for Data, Same event ============

  void processData(TableCollision const& collision,
                   aod::HfCorrSelEl const& electron,
                   TableTracks const& tracks)
  {
    fillCorrelation(collision, electron, tracks);
  }

  PROCESS_SWITCH(HfCorrelatorHfeHadrons, processData, "Process for Data", false);

  // =======  Process starts for McRec, Same event ============

  void processMcRec(McTableCollision const& mcCollision,
                    aod::HfCorrSelEl const& mcElectron,
                    McTableTracks const& mcTracks)
  {
    fillCorrelation(mcCollision, mcElectron, mcTracks);
  }

  PROCESS_SWITCH(HfCorrelatorHfeHadrons, processMcRec, "Process MC Reco mode", false);

  void processMcGen(McGenTableCollision const& mcCollision, aod::McParticles const& mcParticles, aod::HfMcGenSelEl const& electron)
  {

    BinningTypeMcGen corrBinningMcGen{{zBins, multBinsMcGen}, true};
    int poolBin = corrBinningMcGen.getBin(std::make_tuple(mcCollision.posZ(), mcCollision.multMCFT0A()));

    double ptElectron = 0;
    double phiElectron = 0;
    double etaElectron = 0;
    for (const auto& electronMc : electron) {
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
  PROCESS_SWITCH(HfCorrelatorHfeHadrons, processMcGen, "Process MC Gen  mode", true);
  // ====================== Implement Event mixing on Data ===============================

  // ====================== Implement Event mixing on Data ===================================

  void processDataMixedEvent(TableCollisions const& collision, aod::HfCorrSelEl const& electron, TableTracks const& tracks)
  {
    auto tracksTuple = std::make_tuple(electron, tracks);
    Pair<TableCollisions, aod::HfCorrSelEl, TableTracks, BinningType> pair{corrBinning, 5, -1, collision, tracksTuple, &cache};

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
    Pair<McTableCollisions, aod::HfCorrSelEl, McTableTracks, BinningType> pairMcRec{corrBinning, 5, -1, mccollision, tracksTuple, &cache};

    // loop over the rows of the new table
    for (const auto& [c1, tracks1, c2, tracks2] : pairMcRec) {

      fillMixCorrelation(c1, c2, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(HfCorrelatorHfeHadrons, processMcRecMixedEvent, "Process Mixed Event MC Reco mode", false);
  void processMcGenMixedEvent(McGenTableCollisions const& mcCollision, aod::HfMcGenSelEl const& electrons, aod::McParticles const& mcParticles)
  {

    BinningTypeMcGen corrBinningMcGen{{zBins, multBinsMcGen}, true};

    auto tracksTuple = std::make_tuple(electrons, mcParticles);
    Pair<McGenTableCollisions, aod::HfMcGenSelEl, aod::McParticles, BinningTypeMcGen> pairMcGen{corrBinningMcGen, 5, -1, mcCollision, tracksTuple, &cache};

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
