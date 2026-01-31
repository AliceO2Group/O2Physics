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
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/detail/TypeTruncation.h"
#include <Framework/Output.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include <cstdint>
#include <experimental/type_traits> // required for is_detected
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::math_utils::detail;

#define FLOAT_PRECISION 0xFFFFFFF0
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

namespace o2::aod
{
namespace cfmultiplicity
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float); //! Centrality/multiplicity value
} // namespace cfmultiplicity
DECLARE_SOA_TABLE(CFMultiplicities, "AOD", "CFMULTIPLICITY", cfmultiplicity::Multiplicity); //! Transient multiplicity table

using CFMultiplicity = CFMultiplicities::iterator;
} // namespace o2::aod

struct FilterCF {
  Service<o2::framework::O2DatabasePDG> pdg;

  enum TrackSelectionCuts1 : uint8_t {
    kTrackSelected = BIT(0),
    kITS5Clusters = BIT(1),
    kTPCCrossedRows = BIT(2),
    kTPCClusters = BIT(3),
    kchi2perTPC = BIT(4),
    kchi2perITS = BIT(5),
  };

  enum TrackSelectionCuts2 : uint8_t {
    kPIDProton = BIT(1)
  };

  // Configuration
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 7.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPt, float, 0.5f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutMCPt, float, 0.5f, "Minimal pT for particles")
  O2_DEFINE_CONFIGURABLE(cfgCutMCEta, float, 0.8f, "Eta range for particles")
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 1, "Verbosity level (0 = major, 1 = per collision)")
  O2_DEFINE_CONFIGURABLE(cfgTrigger, int, 7, "Trigger choice: (0 = none, 7 = sel7, 8 = sel8, 9 = sel8 + kNoSameBunchPileup + kIsGoodZvtxFT0vsPV, 10 = sel8 before April, 2024, 11 = sel8 for MC, 12 = sel8 with low occupancy cut, 13 = sel8 + kNoSameBunchPileup + kIsGoodITSLayersAll -- for OO/NeNe) ")
  O2_DEFINE_CONFIGURABLE(cfgMinOcc, int, 0, "minimum occupancy selection")
  O2_DEFINE_CONFIGURABLE(cfgMaxOcc, int, 3000, "maximum occupancy selection")
  O2_DEFINE_CONFIGURABLE(cfgCollisionFlags, uint16_t, aod::collision::CollisionFlagsRun2::Run2VertexerTracks, "Request collision flags if non-zero (0 = off, 1 = Run2VertexerTracks)")
  O2_DEFINE_CONFIGURABLE(cfgTransientTables, bool, false, "Output transient tables for collision and track IDs to enable successive filtering tasks")
  O2_DEFINE_CONFIGURABLE(cfgTrackSelection, int, 0, "Type of track selection (0 = Run 2/3 without systematics | 1 = Run 3 with systematics |  2 = Run 3 with proton pid selection)")
  O2_DEFINE_CONFIGURABLE(cfgMinMultiplicity, float, -1, "Minimum multiplicity considered for filtering (if value positive)")
  O2_DEFINE_CONFIGURABLE(cfgMcSpecialPDGs, std::vector<int>, {}, "Special MC PDG codes to include in the MC primary particle output (additional to charged particles). Empty = charged particles only.") // needed for some neutral particles
  O2_DEFINE_CONFIGURABLE(nsigmaCutTPCProton, float, 3, "proton nsigma TPC")
  O2_DEFINE_CONFIGURABLE(nsigmaCutTOFProton, float, 3, "proton nsigma TOF")
  O2_DEFINE_CONFIGURABLE(ITSProtonselection, bool, false, "flag for ITS proton nsigma selection")
  O2_DEFINE_CONFIGURABLE(nsigmaCutITSProton, float, 3, "proton nsigma ITS")
  O2_DEFINE_CONFIGURABLE(dcaxymax, float, 999.f, "maximum dcaxy of tracks")
  O2_DEFINE_CONFIGURABLE(dcazmax, float, 999.f, "maximum dcaz of tracks")
  O2_DEFINE_CONFIGURABLE(itsnclusters, int, 5, "minimum number of ITS clusters for tracks")
  O2_DEFINE_CONFIGURABLE(tpcncrossedrows, int, 80, "minimum number of TPC crossed rows for tracks")
  O2_DEFINE_CONFIGURABLE(tpcnclusters, int, 50, "minimum number of TPC clusters found")
  O2_DEFINE_CONFIGURABLE(chi2pertpccluster, float, 2.5, "maximum Chi2 / cluster for the TPC track segment")
  O2_DEFINE_CONFIGURABLE(chi2peritscluster, float, 36, "maximum Chi2 / cluster for the ITS track segment")
  O2_DEFINE_CONFIGURABLE(cfgEstimatorBitMask, uint16_t, 0, "BitMask for multiplicity estimators to be included in the CFMultSet tables.");

  // Filters and input definitions
  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter collisionVertexTypeFilter = (cfgCollisionFlags == 0) || ((aod::collision::flags & cfgCollisionFlags) == cfgCollisionFlags);

  // TODO how to have this in the second task? For now they are copied
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt);
  Filter trackSelection = (requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true);

  Filter mcCollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;

  OutputObj<TH3F> yields{TH3F("yields", "centrality vs pT vs eta", 100, 0, 100, 40, 0, 20, 100, -2, 2)};
  OutputObj<TH3F> etaphi{TH3F("etaphi", "centrality vs eta vs phi", 100, 0, 100, 100, -2, 2, 200, 0, 2 * M_PI)};

  HistogramRegistry registrytrackQA{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Produces<aod::CFCollisions> outputCollisions;
  Produces<aod::CFTracks> outputTracks;

  Produces<aod::CFCollLabels> outputMcCollisionLabels;
  Produces<aod::CFTrackLabels> outputTrackLabels;

  Produces<aod::CFMcCollisions> outputMcCollisions;
  Produces<aod::CFMcParticles> outputMcParticles;

  Produces<aod::CFCollRefs> outputCollRefs;
  Produces<aod::CFTrackRefs> outputTrackRefs;
  Produces<aod::CFMcParticleRefs> outputMcParticleRefs;

  Produces<aod::CFMultSets> outputMultSets;
  std::vector<float> multiplicities{};

  // persistent caches
  std::vector<bool> mcReconstructedCache;
  std::vector<int> mcParticleLabelsCache;

  void init(InitContext&)
  {
    if (doprocessTrackQA) {
      registrytrackQA.add("zvtx", "Z Vertex position;  posz (cm); Events", HistType::kTH1F, {{100, -12, 12}});
      registrytrackQA.add("eta", "eta distribution;  eta; arb. units", HistType::kTH1F, {{100, -2, 2}});
      registrytrackQA.add("pT", "pT distribution;  #it{p}_{T} (GeV/#it{c}); arb. units", HistType::kTH1F, {{1000, 0, 30}});
      registrytrackQA.add("ptdcaxy", "pT vs DCAxy;  #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", HistType::kTH2F, {{100, 0, 10}, {200, -1, 1}});
      registrytrackQA.add("ptdcaz", "pT vs DCAz; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", HistType::kTH2F, {{100, 0, 10}, {600, -3.0, 3.0}});
      registrytrackQA.add("tpcxrows", "TPC crossed rows; TPC X-rows; Counts", HistType::kTH1F, {{180, 0, 180}});
      registrytrackQA.add("tpcnclst", "TPC found clusters; TPC N_{cls}; Counts", HistType::kTH1F, {{180, 0, 180}});
      registrytrackQA.add("itsnclst", "ITS clusters; ITS N_{cls}; Counts", HistType::kTH1F, {{10, 0, 10}});
      registrytrackQA.add("chi2tpc", "Chi2 per TPC cluster; #chi^{2}/TPC cluster; Counts", HistType::kTH1F, {{100, 0, 10}});
      registrytrackQA.add("chi2its", "Chi2 per ITS cluster; #chi^{2}/ITS cluster; Counts", HistType::kTH1F, {{60, 0, 60}});
    }
  }

  template <typename TCollision>
  bool keepCollision(TCollision& collision)
  {
    bool isMultSelected = false;
    if (collision.multiplicity() >= cfgMinMultiplicity)
      isMultSelected = true;

    if (cfgTrigger == 0) {
      return true;
    } else if (cfgTrigger == 7) {
      return isMultSelected && collision.alias_bit(kINT7) && collision.sel7();
    } else if (cfgTrigger == 8) {
      return isMultSelected && collision.sel8();
    } else if (cfgTrigger == 9) { // relevant only for Pb-Pb
      return isMultSelected && collision.sel8() && collision.selection_bit(aod::evsel::kNoSameBunchPileup) && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && collision.selection_bit(aod::evsel::kIsGoodITSLayersAll);
    } else if (cfgTrigger == 10) { // TVX trigger only (sel8 selection before April, 2024)
      return isMultSelected && collision.selection_bit(aod::evsel::kIsTriggerTVX);
    } else if (cfgTrigger == 11) { // sel8 selection for MC
      return isMultSelected && collision.selection_bit(aod::evsel::kIsTriggerTVX) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder);
    } else if (cfgTrigger == 12) { // relevant only for Pb-Pb with occupancy cuts and rejection of the collisions which have other events nearby
      int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy >= cfgMinOcc && occupancy < cfgMaxOcc)
        return isMultSelected && collision.sel8() && collision.selection_bit(aod::evsel::kNoSameBunchPileup) && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && collision.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) && collision.selection_bit(aod::evsel::kIsGoodITSLayersAll);
      else
        return false;
    } else if (cfgTrigger == 13) { // relevant for pO/OO/NeNe --recommended by Physics Board on 27.01.2026
      return isMultSelected && collision.sel8() && collision.selection_bit(aod::evsel::kNoSameBunchPileup) && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV);
    }
    return false;
  }

  using TrackType = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCPr, aod::pidTOFPr, aod::TracksDCA>>;

  template <typename T>
  bool selectionPIDProton(const T& candidate)
  {
    o2::aod::ITSResponse itsResponse;

    if (ITSProtonselection && candidate.pt() <= 0.6 && !(itsResponse.nSigmaITS<o2::track::PID::Proton>(candidate) > nsigmaCutITSProton)) {
      return false;
    }
    if (ITSProtonselection && candidate.pt() > 0.6 && candidate.pt() <= 0.8 && !(itsResponse.nSigmaITS<o2::track::PID::Proton>(candidate) > nsigmaCutITSProton)) {
      return false;
    }

    if (candidate.hasTOF()) {
      if (candidate.pt() < 0.7 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPCProton) {
        return true;
      }
      if (candidate.p() >= 0.7 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPCProton && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOFProton) {
        return true;
      }
    } else {
      if (std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPCProton) {
        return true;
      }
    }
    return false;
  }

  template <typename T, typename = void>
  struct HasProtonPID : std::false_type {
  };

  template <typename T>
  struct HasProtonPID<T, std::void_t<
                           decltype(std::declval<T>().tpcNSigmaPr()),
                           decltype(std::declval<T>().tofNSigmaPr()),
                           decltype(std::declval<T>().hasTOF())>>
    : std::true_type {
  };

  template <typename TTrack = TrackType>
  uint8_t getTrackType(const TTrack& track)
  {
    if (cfgTrackSelection == 0) {
      if (track.isGlobalTrack()) {
        return 1;
      } else if (track.isGlobalTrackSDD()) {
        return 2;
      }
      return 0;
    } else if (cfgTrackSelection == 1) {
      uint8_t trackType = 0;
      if (track.isGlobalTrack()) {
        trackType |= kTrackSelected;
        if (track.itsNCls() >= itsnclusters) {
          trackType |= kITS5Clusters;
        }
        if (track.tpcNClsCrossedRows() >= tpcncrossedrows) {
          trackType |= kTPCCrossedRows;
        }
        if (track.tpcNClsFound() >= tpcnclusters) {
          trackType |= kTPCClusters;
        }
        if (track.tpcChi2NCl() <= chi2pertpccluster) {
          trackType |= kchi2perTPC;
        }
        if (track.itsChi2NCl() <= chi2peritscluster) {
          trackType |= kchi2perITS;
        }
      }
      return trackType;
    } else if (cfgTrackSelection == 2) {
      uint8_t trackType = 0;
      if constexpr (HasProtonPID<TTrack>::value) {
        if (track.isGlobalTrack() && (track.itsNCls() >= itsnclusters) && (track.tpcNClsCrossedRows() >= tpcncrossedrows) && selectionPIDProton(track)) {
          trackType |= kPIDProton;
        }
      }
      return trackType;
    }

    LOGF(fatal, "Invalid setting for cfgTrackSelection: %d", cfgTrackSelection.value);
    return 0;
  }

  template <class T>
  using HasMultTables = decltype(std::declval<T&>().multNTracksPV());

  /// \brief Templetized process data for a given collision and its associated tracks
  /// \param collision The collision object containing information about the collision
  /// \param tracks The collection of tracks associated with the collision
  template <typename C1, typename T1>
  void processDataT(const C1& collision, const T1& tracks)
  {
    if (cfgVerbosity > 0) {
      LOGF(info, "processData: Tracks for collision: %d | Vertex: %.1f (%d) | INT7: %d | Multiplicity: %.1f", tracks.size(), collision.posZ(), collision.flags(), collision.sel7(), collision.multiplicity());
    }

    if (!keepCollision(collision)) {
      return;
    }

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    outputCollisions(bc.runNumber(), collision.posZ(), collision.multiplicity(), bc.timestamp());

    if constexpr (std::experimental::is_detected<HasMultTables, C1>::value) {
      multiplicities.clear();
      if (cfgEstimatorBitMask & aod::cfmultset::CentFT0C)
        multiplicities.push_back(collision.centFT0C());
      if (cfgEstimatorBitMask & aod::cfmultset::MultFV0A)
        multiplicities.push_back(collision.multFV0A());
      if (cfgEstimatorBitMask & aod::cfmultset::MultNTracksPV)
        multiplicities.push_back(collision.multNTracksPV());
      if (cfgEstimatorBitMask & aod::cfmultset::MultNTracksGlobal)
        multiplicities.push_back(collision.multNTracksGlobal());
      if (cfgEstimatorBitMask & aod::cfmultset::CentFT0M)
        multiplicities.push_back(collision.centFT0M());
      outputMultSets(multiplicities);
    }

    if (cfgTransientTables)
      outputCollRefs(collision.globalIndex());
    for (auto& track : tracks) {
      if ((std::abs(track.dcaXY()) > dcaxymax) || (std::abs(track.dcaZ()) > dcazmax)) {
        continue;
      }

      outputTracks(outputCollisions.lastIndex(), track.pt(), track.eta(), track.phi(), track.sign(), getTrackType(track));
      if (cfgTransientTables)
        outputTrackRefs(collision.globalIndex(), track.globalIndex());

      yields->Fill(collision.multiplicity(), track.pt(), track.eta());
      etaphi->Fill(collision.multiplicity(), track.eta(), track.phi());
    }
  }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CFMultiplicities>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>> const& tracks)
  {
    processDataT(collision, tracks);
  }
  PROCESS_SWITCH(FilterCF, processData, "Process data", true);

  void processDataPid(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CFMultiplicities>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCPr, aod::pidTOFPr, aod::TracksDCA>> const& tracks)
  {
    processDataT(collision, tracks);
  }
  PROCESS_SWITCH(FilterCF, processDataPid, "Process data with PID", false);

  void processDataMults(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CFMultiplicities, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults, aod::FV0Mults, aod::MultsGlobal>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>> const& tracks)
  {
    processDataT(collision, tracks);
  }
  PROCESS_SWITCH(FilterCF, processDataMults, "Process data with multiplicity sets", false);

  void processTrackQA(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CFMultiplicities>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>> const& tracks)
  {
    registrytrackQA.fill(HIST("zvtx"), collision.posZ());
    for (const auto& track : tracks) {
      if (!track.isGlobalTrack()) {
        return; // trackQA for global tracks only
      }
      registrytrackQA.fill(HIST("eta"), track.eta());
      registrytrackQA.fill(HIST("pT"), track.pt());
      registrytrackQA.fill(HIST("ptdcaxy"), track.pt(), track.dcaXY());
      registrytrackQA.fill(HIST("ptdcaz"), track.pt(), track.dcaZ());
      registrytrackQA.fill(HIST("tpcxrows"), track.tpcNClsCrossedRows());
      registrytrackQA.fill(HIST("tpcnclst"), track.tpcNClsFound());
      registrytrackQA.fill(HIST("itsnclst"), track.itsNCls());
      if (track.tpcNClsFound() > 0)
        registrytrackQA.fill(HIST("chi2tpc"), track.tpcChi2NCl());
      if (track.itsNCls() > 0)
        registrytrackQA.fill(HIST("chi2its"), track.itsChi2NCl());
    }
  }
  PROCESS_SWITCH(FilterCF, processTrackQA, "Process track QA", false);

  /// \brief Process MC data for a given set of MC collisions and associated particles and tracks
  /// \param mcCollisions The collection of MC collisions
  /// \param allParticles The collection of all MC particles
  /// \param allCollisions The collection of all collisions, joined with MC collision labels and
  ///                      event selections
  /// \param tracks The collection of tracks, filtered by selection criteria
  /// \param bcs The collection of bunch crossings with timestamps
  template <typename C1, typename T1>
  void processMCT(aod::McCollisions const& mcCollisions, aod::McParticles const& allParticles,
                  C1 const& allCollisions,
                  T1 const& tracks,
                  aod::BCsWithTimestamps const&)
  {
    mcReconstructedCache.reserve(allParticles.size());
    mcParticleLabelsCache.reserve(allParticles.size());
    mcReconstructedCache.clear();
    mcParticleLabelsCache.clear();
    for (int i = 0; i < allParticles.size(); i++) {
      mcReconstructedCache.push_back(false);
      mcParticleLabelsCache.push_back(-1);
    }

    // PASS 1 on collisions: check which particles are kept
    for (auto& collision : allCollisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "processMC:   Tracks for collision %d: %d | Vertex: %.1f (%d) | INT7: %d", collision.globalIndex(), groupedTracks.size(), collision.posZ(), collision.flags(), collision.sel7());
      }

      if (!keepCollision(collision)) {
        continue;
      }

      for (auto& track : groupedTracks) {
        if (track.has_mcParticle()) {
          mcReconstructedCache[track.mcParticleId()] = true;
        }
      }
    }

    for (auto& mcCollision : mcCollisions) {
      auto particles = allParticles.sliceBy(perMcCollision, mcCollision.globalIndex());

      if (cfgVerbosity > 0) {
        LOGF(info, "processMC: Particles for MC collision %d: %d | Vertex: %.1f", mcCollision.globalIndex(), particles.size(), mcCollision.posZ());
      }

      // Store selected MC particles and MC collisions
      int multiplicity = 0;
      for (auto& particle : particles) {
        int8_t sign = 0;
        TParticlePDG* pdgparticle = pdg->GetParticle(particle.pdgCode());
        if (pdgparticle != nullptr) {
          sign = (pdgparticle->Charge() > 0) ? 1.0 : ((pdgparticle->Charge() < 0) ? -1.0 : 0.0);
        }

        bool special = !cfgMcSpecialPDGs->empty() && std::find(cfgMcSpecialPDGs->begin(), cfgMcSpecialPDGs->end(), particle.pdgCode()) != cfgMcSpecialPDGs->end();
        bool primary = particle.isPhysicalPrimary() && sign != 0 && std::abs(particle.eta()) < cfgCutMCEta && particle.pt() > cfgCutMCPt;
        if (primary) {
          multiplicity++;
        }
        if (mcReconstructedCache[particle.globalIndex()] || primary || special) {
          // keep particle

          // use highest bit to flag if it is reconstructed
          uint8_t flags = particle.flags() & ~aod::cfmcparticle::kReconstructed; // clear bit in case of clashes in the future
          if (mcReconstructedCache[particle.globalIndex()]) {
            flags |= aod::cfmcparticle::kReconstructed;
          }

          // NOTE using "outputMcCollisions.lastIndex()+1" here to allow filling of outputMcCollisions *after* the loop
          outputMcParticles(outputMcCollisions.lastIndex() + 1, truncateFloatFraction(particle.pt(), FLOAT_PRECISION), truncateFloatFraction(particle.eta(), FLOAT_PRECISION),
                            truncateFloatFraction(particle.phi(), FLOAT_PRECISION), sign, particle.pdgCode(), flags);
          if (cfgTransientTables)
            outputMcParticleRefs(outputMcCollisions.lastIndex() + 1, particle.globalIndex());

          // relabeling array
          mcParticleLabelsCache[particle.globalIndex()] = outputMcParticles.lastIndex();
        }
      }
      outputMcCollisions(mcCollision.posZ(), multiplicity);
    }

    // PASS 2 on collisions: store collisions and tracks
    for (auto& collision : allCollisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "processMC:   Tracks for collision %d: %d | Vertex: %.1f (%d) | INT7: %d", collision.globalIndex(), groupedTracks.size(), collision.posZ(), collision.flags(), collision.sel7());
      }

      if (!keepCollision(collision)) {
        continue;
      }

      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      // NOTE works only when we store all MC collisions (as we do here)
      outputCollisions(bc.runNumber(), collision.posZ(), collision.multiplicity(), bc.timestamp());
      outputMcCollisionLabels(collision.mcCollisionId());

      if constexpr (std::experimental::is_detected<HasMultTables, C1>::value) {
        multiplicities.clear();
        if (cfgEstimatorBitMask & aod::cfmultset::CentFT0C)
          multiplicities.push_back(collision.centFT0C());
        if (cfgEstimatorBitMask & aod::cfmultset::MultFV0A)
          multiplicities.push_back(collision.multFV0A());
        if (cfgEstimatorBitMask & aod::cfmultset::MultNTracksPV)
          multiplicities.push_back(collision.multNTracksPV());
        if (cfgEstimatorBitMask & aod::cfmultset::MultNTracksGlobal)
          multiplicities.push_back(collision.multNTracksGlobal());
        if (cfgEstimatorBitMask & aod::cfmultset::CentFT0M)
          multiplicities.push_back(collision.centFT0M());
        outputMultSets(multiplicities);
      }

      if (cfgTransientTables)
        outputCollRefs(collision.globalIndex());

      for (auto& track : groupedTracks) {
        int mcParticleId = track.mcParticleId();
        if (mcParticleId >= 0) {
          mcParticleId = mcParticleLabelsCache[track.mcParticleId()];
          if (mcParticleId < 0) {
            LOGP(fatal, "processMC:     Track {} is referring to a MC particle which we do not store {} {} (reco flag {})", track.index(), track.mcParticleId(), mcParticleId, static_cast<bool>(mcReconstructedCache[track.mcParticleId()]));
          }
        }
        outputTracks(outputCollisions.lastIndex(), truncateFloatFraction(track.pt()), truncateFloatFraction(track.eta()), truncateFloatFraction(track.phi()), track.sign(), getTrackType(track));
        outputTrackLabels(mcParticleId);
        if (cfgTransientTables)
          outputTrackRefs(collision.globalIndex(), track.globalIndex());

        yields->Fill(collision.multiplicity(), track.pt(), track.eta());
        etaphi->Fill(collision.multiplicity(), track.eta(), track.phi());
      }
    }
  }

  // NOTE not filtering collisions here because in that case there can be tracks referring to MC particles which are not part of the selected MC collisions
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  void processMC(aod::McCollisions const& mcCollisions, aod::McParticles const& allParticles,
                 soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CFMultiplicities> const& allCollisions,
                 soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels, aod::TrackSelection>> const& tracks,
                 aod::BCsWithTimestamps const& bcs)
  {
    processMCT(mcCollisions, allParticles, allCollisions, tracks, bcs);
  }
  PROCESS_SWITCH(FilterCF, processMC, "Process MC", false);

  // NOTE not filtering collisions here because in that case there can be tracks referring to MC particles which are not part of the selected MC collisions
  void processMCPid(aod::McCollisions const& mcCollisions, aod::McParticles const& allParticles,
                    soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CFMultiplicities> const& allCollisions,
                    soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels, aod::TrackSelection, aod::pidTPCPr, aod::pidTOFPr, aod::TracksDCA>> const& tracks,
                    aod::BCsWithTimestamps const& bcs)
  {
    processMCT(mcCollisions, allParticles, allCollisions, tracks, bcs);
  }
  PROCESS_SWITCH(FilterCF, processMCPid, "Process MC with PID", false);

  void processMCMults(aod::McCollisions const& mcCollisions, aod::McParticles const& allParticles,
                      soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CFMultiplicities, aod::CentFT0Cs, aod::PVMults, aod::FV0Mults, aod::MultsGlobal> const& allCollisions,
                      soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels, aod::TrackSelection>> const& tracks,
                      aod::BCsWithTimestamps const& bcs)
  {
    processMCT(mcCollisions, allParticles, allCollisions, tracks, bcs);
  }

  PROCESS_SWITCH(FilterCF, processMCMults, "Process MC with multiplicity sets", false);

  void processMCGen(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& particles)
  {
    float multiplicity = 0.0f;
    for (auto& particle : particles) {
      if (!particle.isPhysicalPrimary() || std::abs(particle.eta()) > cfgCutMCEta || particle.pt() < cfgCutMCPt)
        continue;
      int8_t sign = 0;
      if (TParticlePDG* pdgparticle = pdg->GetParticle(particle.pdgCode()))
        if ((sign = pdgparticle->Charge()) != 0)
          multiplicity += 1.0f;
      outputMcParticles(outputMcCollisions.lastIndex() + 1, truncateFloatFraction(particle.pt(), FLOAT_PRECISION),
                        truncateFloatFraction(particle.eta(), FLOAT_PRECISION),
                        truncateFloatFraction(particle.phi(), FLOAT_PRECISION),
                        sign, particle.pdgCode(), particle.flags());
    }
    outputMcCollisions(mcCollision.posZ(), multiplicity);
  }
  PROCESS_SWITCH(FilterCF, processMCGen, "Process MCGen", false);
};

struct MultiplicitySelector {
  Produces<aod::CFMultiplicities> output;

  O2_DEFINE_CONFIGURABLE(cfgCutPt, float, 0.5f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")

  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt);
  Filter trackSelection = (requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true);

  void init(InitContext&)
  {
    int enabledFunctions = 0;
    if (doprocessRun2V0M) {
      enabledFunctions++;
    }
    if (doprocessTracks) {
      enabledFunctions++;
    }
    if (doprocessFT0M) {
      enabledFunctions++;
    }
    if (doprocessFT0C) {
      enabledFunctions++;
    }
    if (doprocessFT0CVariant1) {
      enabledFunctions++;
    }
    if (doprocessFT0A) {
      enabledFunctions++;
    }
    if (doprocessCentNGlobal) {
      enabledFunctions++;
    }
    if (doprocessMCGen) {
      enabledFunctions++;
    }

    if (enabledFunctions != 1) {
      LOGP(fatal, "{} multiplicity selectors enabled but we need exactly 1.", enabledFunctions);
    }
  }

  void processTracks(aod::Collision const&, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
  {
    output(tracks.size());
  }
  PROCESS_SWITCH(MultiplicitySelector, processTracks, "Select track count as multiplicity", false);

  void processFT0M(aod::CentFT0Ms const& centralities)
  {
    for (auto& c : centralities) {
      output(c.centFT0M());
    }
  }
  PROCESS_SWITCH(MultiplicitySelector, processFT0M, "Select FT0M centrality as multiplicity", false);

  void processFT0C(aod::CentFT0Cs const& centralities)
  {
    for (auto& c : centralities) {
      output(c.centFT0C());
    }
  }
  PROCESS_SWITCH(MultiplicitySelector, processFT0C, "Select FT0C centrality as multiplicity", false);

  void processFT0CVariant1(aod::CentFT0CVariant1s const& centralities)
  {
    for (auto& c : centralities) {
      output(c.centFT0CVariant1());
    }
  }
  PROCESS_SWITCH(MultiplicitySelector, processFT0CVariant1, "Select FT0CVariant1 centrality as multiplicity", false);

  void processFT0A(aod::CentFT0As const& centralities)
  {
    for (auto& c : centralities) {
      output(c.centFT0A());
    }
  }
  PROCESS_SWITCH(MultiplicitySelector, processFT0A, "Select FT0A centrality as multiplicity", false);

  void processCentNGlobal(aod::CentNGlobals const& centralities)
  {
    for (auto& c : centralities) {
      output(c.centNGlobal());
    }
  }
  PROCESS_SWITCH(MultiplicitySelector, processCentNGlobal, "Select CentNGlobal centrality as multiplicity", false);

  void processRun2V0M(aod::CentRun2V0Ms const& centralities)
  {
    for (auto& c : centralities) {
      output(c.centRun2V0M());
    }
  }
  PROCESS_SWITCH(MultiplicitySelector, processRun2V0M, "Select V0M centrality as multiplicity", true);

  void processMCGen(aod::McCollision const&, aod::McParticles const& particles)
  {
    output(particles.size());
  }
  PROCESS_SWITCH(MultiplicitySelector, processMCGen, "Select MC particle count as multiplicity", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FilterCF>(cfgc),
    adaptAnalysisTask<MultiplicitySelector>(cfgc)};
}
