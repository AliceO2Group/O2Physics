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

/// \file nucleibalance.cxx
/// \brief task for the balance function and correlations for nuclei for O2 analysis. First part is inspired from PWGCF/Tasks/correlations.cxx
/// \author Sushanta Tripathy <sushanta.tripathy@cern.ch>

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TFormula.h>
#include <TH1F.h>
#include <THn.h>
#include <THnSparse.h>
#include <TPDGCode.h>
#include <TVector2.h>

#include <cmath>
#include <cstring>
#include <deque>
#include <experimental/type_traits>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

static constexpr float PairCutOff = -1.f;
static constexpr float CfgPairCutDefaults[1][5] = {{PairCutOff, PairCutOff, PairCutOff, PairCutOff, PairCutOff}};

struct Nucleibalance {
  SliceCache cache;

  // Configuration
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 7.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPt, float, 0.5f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutMCPt, float, 0.5f, "Minimal pT for MC particles (AO2D-MC mode)")
  O2_DEFINE_CONFIGURABLE(cfgCutMCEta, float, 0.8f, "Eta range for MC particles (AO2D-MC mode)")
  O2_DEFINE_CONFIGURABLE(cfgUseFT0M, int, 1, "Use FT0M centrality (0-100) as multiplicity axis (1=ON, 0=use Ntracks/multiplicity())")
  // Track-quality options (AO2D mode). Default selection corresponds to global tracks.
  O2_DEFINE_CONFIGURABLE(cfgTPCNClsMin, int, 70, "Minimum number of TPC clusters (tpcNClsFound) in AO2D mode")
  O2_DEFINE_CONFIGURABLE(cfgDcaXYMax, float, 0.1f, "Max |DCA_{xy}| to PV (cm) in AO2D mode")
  O2_DEFINE_CONFIGURABLE(cfgDcaZMax, float, 0.2f, "Max |DCA_{z}| to PV (cm) in AO2D mode")
  O2_DEFINE_CONFIGURABLE(chi2pertpccluster, float, 2.5f, "Maximum Chi2/cluster for the TPC track segment in AO2D mode")
  O2_DEFINE_CONFIGURABLE(chi2peritscluster, float, 36.f, "Maximum Chi2/cluster for the ITS track segment in AO2D mode")
  O2_DEFINE_CONFIGURABLE(itsnclusters, int, 5, "Minimum number of ITS clusters in AO2D mode")

  O2_DEFINE_CONFIGURABLE(cfgPtOrder, int, 1, "Only consider pairs for which pT,1 < pT,2 (0 = OFF, 1 = ON)");
  O2_DEFINE_CONFIGURABLE(cfgTriggerCharge, int, 0, "Select on charge of trigger particle: 0 = all; 1 = positive; -1 = negative");
  O2_DEFINE_CONFIGURABLE(cfgAssociatedCharge, int, 0, "Select on charge of associated particle: 0 = all charged; 1 = positive; -1 = negative");
  O2_DEFINE_CONFIGURABLE(cfgPairCharge, int, 0, "Select on charge of particle pair: 0 = all; 1 = like sign; -1 = unlike sign");

  O2_DEFINE_CONFIGURABLE(cfgTwoTrackCut, float, -1, "Two track cut: -1 = off; >0 otherwise distance value (suggested: 0.02)");
  O2_DEFINE_CONFIGURABLE(cfgTwoTrackCutMinRadius, float, 0.8f, "Two track cut: radius in m from which two track cuts are applied");
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, int, 0, "0 = OFF and 1 = ON for local efficiency");
  O2_DEFINE_CONFIGURABLE(cfgCentBinsForMC, int, 0, "0 = OFF and 1 = ON for data like multiplicity/centrality bins for MC steps");
  O2_DEFINE_CONFIGURABLE(cfgTrackBitMask, uint16_t, 1, "BitMask for track selection systematics; refer to the enum TrackSelectionCuts in filtering task (default=1 selects global tracks)");
  O2_DEFINE_CONFIGURABLE(cfgMultCorrelationsMask, uint16_t, 0, "Selection bitmask for the multiplicity correlations. This should match the filter selection cfgEstimatorBitMask.")
  O2_DEFINE_CONFIGURABLE(cfgMultCutFormula, std::string, "", "Multiplicity correlations cut formula. A result greater than zero results in accepted event. Parameters: [cFT0C] FT0C centrality, [mFV0A] V0A multiplicity, [mGlob] global track multiplicity, [mPV] PV track multiplicity")

  // PID and species selection for AO2D-based correlations (pi, K, p, d)
  O2_DEFINE_CONFIGURABLE(cfgUseTPCOnlyPID, int, 1, "Use only TPC PID (1 = TPC only, 0 = require both TPC and TOF when available)");
  O2_DEFINE_CONFIGURABLE(cfgNsigmaTPCPi, float, 3.0f, "|n#sigma^{TPC}_{#pi}| cut");
  O2_DEFINE_CONFIGURABLE(cfgNsigmaTPCKa, float, 3.0f, "|n#sigma^{TPC}_{K}| cut");
  O2_DEFINE_CONFIGURABLE(cfgNsigmaTPCPr, float, 3.0f, "|n#sigma^{TPC}_{p}| cut");
  O2_DEFINE_CONFIGURABLE(cfgNsigmaTPCDe, float, 3.0f, "|n#sigma^{TPC}_{d}| cut");
  O2_DEFINE_CONFIGURABLE(cfgNsigmaTOFPi, float, 3.0f, "|n#sigma^{TOF}_{#pi}| cut");
  O2_DEFINE_CONFIGURABLE(cfgNsigmaTOFKa, float, 3.0f, "|n#sigma^{TOF}_{K}| cut");
  O2_DEFINE_CONFIGURABLE(cfgNsigmaTOFPr, float, 3.0f, "|n#sigma^{TOF}_{p}| cut");
  O2_DEFINE_CONFIGURABLE(cfgNsigmaTOFDe, float, 3.0f, "|n#sigma^{TOF}_{d}| cut");

  // Species choice for trigger/associated in the BF:
  // 0 = pion, 1 = kaon, 2 = proton, 3 = deuteron, -1 = all charged tracks
  O2_DEFINE_CONFIGURABLE(cfgTriggerSpecies, int, 3, "Trigger species for BF: 0 = #pi, 1 = K, 2 = p, 3 = d, -1 = all charged tracks");
  O2_DEFINE_CONFIGURABLE(cfgAssociatedSpecies, int, 2, "Associated species for BF: 0 = #pi, 1 = K, 2 = p, 3 = d, -1 = all charged tracks");

  // Suggested values: Photon: 0.004; K0 and Lambda: 0.005
  Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut", {CfgPairCutDefaults[0], 5, {"Photon", "K0", "Lambda", "Phi", "Rho"}}, "Pair cuts on various particles"};

  O2_DEFINE_CONFIGURABLE(cfgEfficiencyTrigger, std::string, "", "CCDB path to efficiency object for trigger particles")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyAssociated, std::string, "", "CCDB path to efficiency object for associated particles")

  O2_DEFINE_CONFIGURABLE(cfgNoMixedEvents, int, 5, "Number of mixed events per event")

  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 1, "Verbosity level (0 = major, 1 = per collision)")

  O2_DEFINE_CONFIGURABLE(cfgMcTriggerPDGs, std::vector<int>, {}, "MC PDG codes to use exclusively as trigger particles and exclude from associated particles. Empty = no selection.")

  ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  Configurable<int> cfgTrigger{"cfgTrigger", 4, "Event trigger selection: 0=none, 1=sel8, 2=sel8+(kNoSameBunchPileup+kIsGoodZvtxFT0vsPV+kIsGoodITSLayersAll), 3=sel8+occupancy+(kNoCollInTimeRangeStandard)+(kNoSameBunchPileup+kIsGoodZvtxFT0vsPV+kIsGoodITSLayersAll), 4=sel8+(kNoSameBunchPileup+kIsGoodZvtxFT0vsPV)"};
  Configurable<int> cfgMinOcc{"cfgMinOcc", 0, "minimum occupancy selection (for cfgTrigger==3)"};
  Configurable<int> cfgMaxOcc{"cfgMaxOcc", 3000, "maximum occupancy selection (for cfgTrigger==3)"};

  // Named trigger codes to avoid magic numbers in keepCollisionAO2D
  // ---- Ion/nucleus PDG encoding helpers (10LZZZAAAI) ----
  // Note: these are *format* constants (not particle PDG species codes)
  static constexpr int IonCodeThreshold = 1000000000; // 10^9
  static constexpr int IonZDivisor = 10000;
  static constexpr int IonZModulo = 1000;
  static constexpr int PdgElectron = static_cast<int>(PDG_t::kElectron);
  static constexpr int PdgMuon = static_cast<int>(PDG_t::kMuonMinus);
  static constexpr int PdgPion = static_cast<int>(PDG_t::kPiPlus);
  static constexpr int PdgKaon = static_cast<int>(PDG_t::kKPlus);
  static constexpr int PdgProton = static_cast<int>(PDG_t::kProton);
  static constexpr int TriggerNone = 0;
  static constexpr int TriggerSel8 = 1;
  static constexpr int TriggerSel8Quality = 2;
  static constexpr int TriggerSel8OccQuality = 3;
  static constexpr int TriggerSel8NoSbpZvtx = 4;
  template <typename TCollision>
  bool keepCollisionAO2D(TCollision const& collision) const
  {
    if (cfgTrigger.value == TriggerNone) {
      return true;
    } else if (cfgTrigger.value == TriggerSel8) {
      return collision.sel8();
    } else if (cfgTrigger.value == TriggerSel8Quality) {
      return collision.sel8() &&
             collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
             collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
             collision.selection_bit(aod::evsel::kIsGoodITSLayersAll);
    } else if (cfgTrigger.value == TriggerSel8OccQuality) {
      const int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < cfgMinOcc.value || occupancy >= cfgMaxOcc.value) {
        return false;
      }
      return collision.sel8() &&
             collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
             collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
             collision.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) &&
             collision.selection_bit(aod::evsel::kIsGoodITSLayersAll);
    } else if (cfgTrigger.value == TriggerSel8NoSbpZvtx) {
      return collision.sel8() &&
             collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
             collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV);
    }

    LOGF(warn, "Invalid cfgTrigger=%d. Accepting all collisions.", cfgTrigger.value);
    return true;
  }

  template <typename TCollision>
  float eventMultiplicityFT0MOrFallback(TCollision const& collision, float fallback) const
  {
    if (cfgUseFT0M.value == 0) {
      return fallback;
    }

    // Prefer FT0M centrality if present in this collision table
    if constexpr (requires { collision.centFT0M(); }) {
      const float v = collision.centFT0M();
      if (v >= 0.f) {
        return v; // expected 0..100
      }
    }

    // Some tables may expose a validity bit
    if constexpr (requires { collision.centFT0MValid(); }) {
      if (collision.centFT0MValid()) {
        if constexpr (requires { collision.centFT0M(); }) {
          return collision.centFT0M();
        }
      }
    }

    return fallback;
  }

  static int chargeFromPdg(int pdg)
  {
    const int apdg = std::abs(pdg);

    // Ions/nuclei: PDG code format 10LZZZAAAI -> Z is encoded in digits [7..5]
    if (apdg >= IonCodeThreshold) {
      const int z = (apdg / IonZDivisor) % IonZModulo;
      return (pdg >= 0) ? z : -z;
    }

    // Common charged hadrons/leptons (extend if needed)
    switch (apdg) {
      case PdgElectron: // e
      case PdgMuon:     // mu
      case PdgPion:     // pi
      case PdgKaon:     // K
      case PdgProton:   // p
        return (pdg >= 0) ? 1 : -1;
      default:
        return 0;
    }
  }

  // This filter is applied to AOD and derived data (column names are identical)
  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  // This filter is only applied to AOD
  Filter collisionVertexTypeFilter = (aod::collision::flags & static_cast<uint16_t>(aod::collision::CollisionFlagsRun2::Run2VertexerTracks)) == static_cast<uint16_t>(aod::collision::CollisionFlagsRun2::Run2VertexerTracks);

  // Track filters
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));
  Filter cfTrackFilter = (nabs(aod::cftrack::eta) < cfgCutEta) && (aod::cftrack::pt > cfgCutPt) && ((aod::track::trackType & (uint8_t)cfgTrackBitMask) == (uint8_t)cfgTrackBitMask);

  // MC filters
  Filter cfMCCollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  Filter cfMCParticleFilter = (nabs(aod::cfmcparticle::eta) < cfgCutEta) && (aod::cfmcparticle::pt > cfgCutPt); // && (aod::cfmcparticle::sign != 0); //check the sign manually, some specials may be neutral

  // Output definitions
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  // persistent caches
  std::vector<float> efficiencyAssociatedCache;

  std::unique_ptr<TFormula> multCutFormula;
  std::array<uint, 4> multCutFormulaParamIndex;

  struct Config {
    bool mPairCuts = false;
    THn* mEfficiencyTrigger = nullptr;
    THn* mEfficiencyAssociated = nullptr;
    bool efficiencyLoaded = false;
  } cfg;

  HistogramRegistry registry{"registry"};
  PairCuts mPairCuts;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // AO2D-based tracks with PID for pi / K / p / d
  using TracksPID = soa::Join<aod::Tracks,
                              aod::TracksExtra,
                              aod::TrackSelection,
                              aod::TracksDCA,
                              aod::pidTPCFullPi,
                              aod::pidTPCFullKa,
                              aod::pidTPCFullPr,
                              aod::pidTPCFullDe,
                              aod::pidTOFFullPi,
                              aod::pidTOFFullKa,
                              aod::pidTOFFullPr,
                              aod::pidTOFFullDe>;
  using TracksPIDFiltered = soa::Filtered<TracksPID>;

  using TracksPIDMC = soa::Join<aod::Tracks,
                                aod::TracksExtra,
                                aod::TrackSelection,
                                aod::TracksDCA,
                                aod::pidTPCFullPi,
                                aod::pidTPCFullKa,
                                aod::pidTPCFullPr,
                                aod::pidTPCFullDe,
                                aod::pidTOFFullPi,
                                aod::pidTOFFullKa,
                                aod::pidTOFFullPr,
                                aod::pidTOFFullDe,
                                aod::McTrackLabels>;

  using CollisionsAO2DMC = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::McCollisionLabels>;

  // group MC particles by MC collision
  Preslice<aod::McParticles> mcParticlesPerCollision = aod::mcparticle::mcCollisionId;

  // Helper: access a row of a soa::Join by row index (operator[] is not available for joins)
  template <typename TTracks>
  static auto trackAt(TTracks const& tracks, uint32_t idx)
  {
    if constexpr (requires { tracks.iteratorAt(idx); }) {
      return tracks.iteratorAt(idx);
    } else {
      return tracks.rawIteratorAt(idx);
    }
  }

  // Helper: AO2D track-quality selection (default: global tracks)
  template <typename TTrack>
  bool passTrackQualityAO2D(const TTrack& trk) const
  {
    // Default: require global tracks when the column exists
    if constexpr (requires { trk.isGlobalTrack(); }) {
      if (!trk.isGlobalTrack()) {
        return false;
      }
    } else if constexpr (requires { trk.isGlobalTrackSDD(); }) {
      // fallback (older tables)
      if (!trk.isGlobalTrackSDD()) {
        return false;
      }
    }

    if constexpr (requires { trk.itsNCls(); }) {
      if (itsnclusters.value > 0 && trk.itsNCls() < itsnclusters.value) {
        return false;
      }
    }

    if constexpr (requires { trk.tpcNClsFound(); }) {
      if (cfgTPCNClsMin.value > 0 && trk.tpcNClsFound() < cfgTPCNClsMin.value) {
        return false;
      }
    }

    if constexpr (requires { trk.tpcChi2NCl(); }) {
      if (chi2pertpccluster.value > 0.f && trk.tpcChi2NCl() > chi2pertpccluster.value) {
        return false;
      }
    }

    if constexpr (requires { trk.itsChi2NCl(); }) {
      if (chi2peritscluster.value > 0.f && trk.itsChi2NCl() > chi2peritscluster.value) {
        return false;
      }
    }

    if constexpr (requires { trk.dcaXY(); }) {
      if (cfgDcaXYMax.value > 0.f && std::abs(trk.dcaXY()) > cfgDcaXYMax.value) {
        return false;
      }
    }

    if constexpr (requires { trk.dcaZ(); }) {
      if (cfgDcaZMax.value > 0.f && std::abs(trk.dcaZ()) > cfgDcaZMax.value) {
        return false;
      }
    }

    return true;
  }

  struct SimpleTrack {
    float eta;
    float phi;
    float pt;
    int charge;
  };

  struct MixEventEntry {
    float multiplicity;
    float zvtx;
    std::vector<SimpleTrack> triggerTracks;
    std::vector<SimpleTrack> associatedTracks;
  };

  // Very simple mixing buffer: keep last cfgNoMixedEvents events
  std::deque<MixEventEntry> mMixEvents;
  std::deque<MixEventEntry> mMixEventsMC;

  // Preslice to group AO2D tracks by collision

  using DerivedCollisions = soa::Filtered<aod::CFCollisions>;
  using DerivedTracks = soa::Filtered<aod::CFTracks>;
  using DerivedTracksWithRefs = soa::Filtered<soa::Join<aod::CFTracks, aod::CFTrackRefs>>;

  void init(o2::framework::InitContext&)
  {
    // --- HISTOGRAMS ---
    registry.add("yields", "multiplicity/centrality vs pT vs eta", {HistType::kTH3F, {{100, 0, 100, "/multiplicity/centrality"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("etaphi", "multiplicity/centrality vs eta vs phi", {HistType::kTH3F, {{100, 0, 100, "multiplicity/centrality"}, {100, -2, 2, "#eta"}, {200, 0, o2::constants::math::TwoPI, "#varphi"}}});

    if (doprocessSameDerivedMultSet) {
      if (cfgMultCorrelationsMask == 0)
        LOGF(fatal, "cfgMultCorrelationsMask can not be 0 when MultSet process functions are in use.");
      std::vector<AxisSpec> multAxes;
      if (cfgMultCorrelationsMask & aod::cfmultset::CentFT0C)
        multAxes.emplace_back(100, 0, 100, "FT0C centrality");
      if (cfgMultCorrelationsMask & aod::cfmultset::MultFV0A)
        multAxes.emplace_back(1000, 0, 100000, "V0A multiplicity");
      if (cfgMultCorrelationsMask & aod::cfmultset::MultNTracksPV)
        multAxes.emplace_back(100, 0, 1000, "Nch PV");
      if (cfgMultCorrelationsMask & aod::cfmultset::MultNTracksGlobal)
        multAxes.emplace_back(100, 0, 1000, "Nch Global");
      registry.add("multCorrelations", "Multiplicity correlations", {HistType::kTHnSparseF, multAxes});
    }
    registry.add("multiplicity", "event multiplicity", {HistType::kTH1F, {{1000, 0, 100, "/multiplicity/centrality"}}});
    registry.add("yvspt", "y vs pT", {HistType::kTH2F, {{100, -1, 1, "y"}, {100, 0, 20, "p_{T}"}}}); // y vs pT for all tracks (control histogram)

    const int maxMixBin = AxisSpec(axisMultiplicity).getNbins() * AxisSpec(axisVertex).getNbins();
    // The bin numbers for the control histograms (eventcount_*) come from getBin(...) and are the following: #mult_bin * #number_of_z_bins + #zbin
    registry.add("eventcount_same", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("eventcount_mixed", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("trackcount_same", "bin", {HistType::kTH2F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {10, -0.5, 9.5}}});
    registry.add("trackcount_mixed", "bin", {HistType::kTH3F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {10, -0.5, 9.5}, {10, -0.5, 9.5}}});

    mPairCuts.SetHistogramRegistry(&registry);

    if (cfgPairCut->get("Photon") > 0 || cfgPairCut->get("K0") > 0 || cfgPairCut->get("Lambda") > 0 || cfgPairCut->get("Phi") > 0 || cfgPairCut->get("Rho") > 0) {
      mPairCuts.SetPairCut(PairCuts::Photon, cfgPairCut->get("Photon"));
      mPairCuts.SetPairCut(PairCuts::K0, cfgPairCut->get("K0"));
      mPairCuts.SetPairCut(PairCuts::Lambda, cfgPairCut->get("Lambda"));
      mPairCuts.SetPairCut(PairCuts::Phi, cfgPairCut->get("Phi"));
      mPairCuts.SetPairCut(PairCuts::Rho, cfgPairCut->get("Rho"));
      cfg.mPairCuts = true;
    }

    if (cfgTwoTrackCut > 0) {
      mPairCuts.SetTwoTrackCuts(cfgTwoTrackCut, cfgTwoTrackCutMinRadius);
    }

    // --- OBJECT INIT ---

    if (!cfgMultCutFormula.value.empty()) {
      multCutFormula = std::make_unique<TFormula>("multCutFormula", cfgMultCutFormula.value.c_str());
      std::fill_n(multCutFormulaParamIndex.begin(), std::size(multCutFormulaParamIndex), ~0u);
      std::array<std::string, 4> pars = {"cFT0C", "mFV0A", "mPV", "mGlob"}; // must correspond the order of MultiplicityEstimators
      for (uint i = 0, n = multCutFormula->GetNpar(); i < n; ++i) {
        auto m = std::find(pars.begin(), pars.end(), multCutFormula->GetParName(i));
        if (m == pars.end()) {

          LOGF(warning, "Unknown parameter in cfgMultCutFormula: %s", multCutFormula->GetParName(i));
          continue;
        }
        if ((cfgMultCorrelationsMask.value & (1u << i)) == 0) {
          LOGF(warning, "The centrality/multiplicity estimator %s is not available to be used in cfgMultCutFormula. Ensure cfgMultCorrelationsMask is correct and matches the CFMultSets in derived data.");
        } else {
          multCutFormulaParamIndex[std::distance(pars.begin(), m)] = i;
          LOGF(info, "Multiplicity cut parameter %s in use.", m->c_str());
        }
      }
    }

    std::vector<AxisSpec> corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity / centrality"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    std::vector<AxisSpec> userAxis;
    std::vector<AxisSpec> userMixingAxis;

    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userMixingAxis));

    same->setTrackEtaCut(cfgCutEta);
    mixed->setTrackEtaCut(cfgCutEta);

    if (!cfgEfficiencyAssociated.value.empty())
      efficiencyAssociatedCache.reserve(512);

    // o2-ccdb-upload -p Users/jgrosseo/correlations/LHC15o -f /tmp/correction_2011_global.root -k correction

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    // static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      // grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <class T>
  using HasMultSet = decltype(std::declval<T&>().multiplicities());

  template <typename TCollision, typename TTracks>
  void fillQA(const TCollision& /*collision*/, float multiplicity, const TTracks& tracks)
  {
    registry.fill(HIST("multiplicity"), multiplicity);
    for (const auto& track1 : tracks) {
      registry.fill(HIST("yields"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("etaphi"), multiplicity, track1.eta(), track1.phi());
    }
  }

  template <class T>
  using HasInvMass = decltype(std::declval<T&>().invMass());
  template <class T>
  using HasPDGCode = decltype(std::declval<T&>().pdgCode());

  template <typename TCollision, typename TTracks1, typename TTracks2>
  void fillQA(const TCollision& collision, float multiplicity, const TTracks1& tracks1, const TTracks2& tracks2)
  {
    for (const auto& track1 : tracks1) {
      registry.fill(HIST("yieldsTrigger"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("etaphiTrigger"), multiplicity, track1.eta(), track1.phi());
    }
    fillQA(collision, multiplicity, tracks2);
  }

  template <CorrelationContainer::CFStep step, typename TTrack>
  bool checkObject(TTrack& track)
  {
    if constexpr (step <= CorrelationContainer::kCFStepAnaTopology) {
      return track.isPhysicalPrimary();
    } else if constexpr (step == CorrelationContainer::kCFStepTrackedOnlyPrim) {
      return track.isPhysicalPrimary() && (track.flags() & aod::cfmcparticle::kReconstructed);
    } else if constexpr (step == CorrelationContainer::kCFStepTracked) {
      return (track.flags() & aod::cfmcparticle::kReconstructed);
    }

    return true;
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracks1, typename TTracks2>
  void fillCorrelations(TTarget target, TTracks1& tracks1, TTracks2& tracks2, float multiplicity, float posZ, int magField, float eventWeight)
  {
    // Helper lambda for pair charge selection
    auto passPairCharge = [this](auto const& t1, auto const& t2) {
      if (cfgPairCharge.value == 0) {
        return true;
      }

      int q1 = 0;
      int q2 = 0;

      if constexpr (requires { t1.sign(); }) {
        q1 = t1.sign();
      } else if constexpr (requires { t1.charge(); }) {
        q1 = t1.charge();
      } else if constexpr (requires { t1.pdgCode(); }) {
        q1 = chargeFromPdg(t1.pdgCode());
      }

      if constexpr (requires { t2.sign(); }) {
        q2 = t2.sign();
      } else if constexpr (requires { t2.charge(); }) {
        q2 = t2.charge();
      } else if constexpr (requires { t2.pdgCode(); }) {
        q2 = chargeFromPdg(t2.pdgCode());
      }

      if (q1 == 0 || q2 == 0) {
        // If we cannot determine both charges, reject the pair for pair-charge selections
        return false;
      }

      const int pairSign = q1 * q2;
      if (cfgPairCharge.value == 1) { // like-sign pairs only
        return pairSign > 0;
      }
      if (cfgPairCharge.value == -1) { // unlike-sign pairs only
        return pairSign < 0;
      }

      return true;
    };

    // Cache efficiency for particles (too many FindBin lookups)
    if constexpr (step == CorrelationContainer::kCFStepCorrected) {
      if (cfg.mEfficiencyAssociated) {
        efficiencyAssociatedCache.clear();
        efficiencyAssociatedCache.reserve(tracks2.size());
        for (const auto& track : tracks2) {
          efficiencyAssociatedCache.push_back(getEfficiencyCorrection(cfg.mEfficiencyAssociated, track.eta(), track.pt(), multiplicity, posZ));
        }
      }
    }

    for (const auto& track1 : tracks1) {
      // LOGF(info, "Track %f | %f | %f  %d %d", track1.eta(), track1.phi(), track1.pt(), track1.isGlobalTrack(), track1.isGlobalTrackSDD());

      if constexpr (step <= CorrelationContainer::kCFStepTracked) {
        if (!checkObject<step>(track1)) {
          continue;
        }
      }

      float triggerWeight = eventWeight;
      if constexpr (step == CorrelationContainer::kCFStepCorrected) {
        if (cfg.mEfficiencyTrigger) {
          triggerWeight *= getEfficiencyCorrection(cfg.mEfficiencyTrigger, track1.eta(), track1.pt(), multiplicity, posZ);
        }
      }

      target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, triggerWeight);

      for (const auto& track2 : tracks2) {
        if constexpr (std::is_same<TTracks1, TTracks2>::value) {
          if (track1.globalIndex() == track2.globalIndex()) {
            // LOGF(info, "Track identical: %f | %f | %f || %f | %f | %f", track1.eta(), track1.phi(), track1.pt(),  track2.eta(), track2.phi(), track2.pt());
            continue;
          }
        }
        if constexpr (step <= CorrelationContainer::kCFStepTracked) {
          if (!checkObject<step>(track2)) {
            continue;
          }
        }

        // Pair charge selection
        if (!passPairCharge(track1, track2)) {
          continue;
        }

        if (cfgPtOrder != 0 && track2.pt() >= track1.pt()) {
          continue;
        }

        if constexpr (std::is_same<TTracks1, TTracks2>::value) {
          if constexpr (step >= CorrelationContainer::kCFStepReconstructed) {
            if (cfg.mPairCuts && mPairCuts.conversionCuts(track1, track2)) {
              continue;
            }
            if (cfgTwoTrackCut > 0 && mPairCuts.twoTrackCut(track1, track2, magField)) {
              continue;
            }
          }
        }

        float associatedWeight = triggerWeight;
        if constexpr (step == CorrelationContainer::kCFStepCorrected) {
          if (cfg.mEfficiencyAssociated) {
            associatedWeight *= efficiencyAssociatedCache[track2.filteredIndex()];
          }
        }

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -o2::constants::math::PIHalf);
        target->getPairHist()->Fill(step, track1.eta() - track2.eta(), track2.pt(), track1.pt(), multiplicity, deltaPhi, posZ, associatedWeight);
      }
    }
  }

  void loadEfficiency(uint64_t timestamp)
  {
    if (cfg.efficiencyLoaded) {
      return;
    }
    if (cfgEfficiencyTrigger.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiencyTrigger.value.c_str(), "READ");
        cfg.mEfficiencyTrigger = reinterpret_cast<THn*>(fEfficiencyTrigger->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyTrigger = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyTrigger, timestamp);
      }
      if (cfg.mEfficiencyTrigger == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiencyTrigger.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for trigger particles from %s (%p)", cfgEfficiencyTrigger.value.c_str(), (void*)cfg.mEfficiencyTrigger);
    }
    if (cfgEfficiencyAssociated.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyAssociated = TFile::Open(cfgEfficiencyAssociated.value.c_str(), "READ");
        cfg.mEfficiencyAssociated = reinterpret_cast<THn*>(fEfficiencyAssociated->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyAssociated = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyAssociated, timestamp);
      }
      if (cfg.mEfficiencyAssociated == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", cfgEfficiencyAssociated.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for associated particles from %s (%p)", cfgEfficiencyAssociated.value.c_str(), (void*)cfg.mEfficiencyAssociated);
    }
    cfg.efficiencyLoaded = true;
  }

  double getEfficiencyCorrection(THn* eff, float eta, float pt, float multiplicity, float posZ)
  {
    int effVars[4];
    effVars[0] = eff->GetAxis(0)->FindBin(eta);
    effVars[1] = eff->GetAxis(1)->FindBin(pt);
    effVars[2] = eff->GetAxis(2)->FindBin(multiplicity);
    effVars[3] = eff->GetAxis(3)->FindBin(posZ);
    return eff->GetBinContent(effVars);
  }

  template <typename TTrack>
  bool passPIDForSpecies(const TTrack& trk, int species)
  {
    // -1 means "all charged tracks"
    if (species < 0) {
      return trk.sign() != 0;
    }

    const bool useTPCOnly = (cfgUseTPCOnlyPID.value != 0);

    switch (species) {
      case 0: { // pion
        if (cfgNsigmaTPCPi.value > 0.f && std::abs(trk.tpcNSigmaPi()) > cfgNsigmaTPCPi.value) {
          return false;
        }
        if (!useTPCOnly && cfgNsigmaTOFPi.value > 0.f && std::abs(trk.tofNSigmaPi()) > cfgNsigmaTOFPi.value) {
          return false;
        }
        return trk.sign() != 0;
      }
      case 1: { // kaon
        if (cfgNsigmaTPCKa.value > 0.f && std::abs(trk.tpcNSigmaKa()) > cfgNsigmaTPCKa.value) {
          return false;
        }
        if (!useTPCOnly && cfgNsigmaTOFKa.value > 0.f && std::abs(trk.tofNSigmaKa()) > cfgNsigmaTOFKa.value) {
          return false;
        }
        return trk.sign() != 0;
      }
      case 2: { // proton
        if (cfgNsigmaTPCPr.value > 0.f && std::abs(trk.tpcNSigmaPr()) > cfgNsigmaTPCPr.value) {
          return false;
        }
        if (!useTPCOnly && cfgNsigmaTOFPr.value > 0.f && std::abs(trk.tofNSigmaPr()) > cfgNsigmaTOFPr.value) {
          return false;
        }
        return trk.sign() != 0;
      }
      case 3: { // deuteron
        if (cfgNsigmaTPCDe.value > 0.f && std::abs(trk.tpcNSigmaDe()) > cfgNsigmaTPCDe.value) {
          return false;
        }
        if (!useTPCOnly && cfgNsigmaTOFDe.value > 0.f && std::abs(trk.tofNSigmaDe()) > cfgNsigmaTOFDe.value) {
          return false;
        }
        return trk.sign() != 0;
      }
      default:
        return false;
    }
  }

  // Simple correlation filler writing directly into CorrelationContainer
  void fillCorrelationsSimple(OutputObj<CorrelationContainer>& target,
                              CorrelationContainer::CFStep step,
                              const std::vector<SimpleTrack>& triggers,
                              const std::vector<SimpleTrack>& associates,
                              float multiplicity,
                              float posZ,
                              float eventWeight)
  {
    auto* trigHist = target->getTriggerHist();
    auto* pairHist = target->getPairHist();

    for (auto const& t : triggers) {
      trigHist->Fill(step, t.pt, multiplicity, posZ, eventWeight);

      for (auto const& a : associates) {
        if (cfgPtOrder != 0 && a.pt >= t.pt) {
          continue;
        }
        // Pair charge selection
        if (cfgPairCharge.value != 0) {
          const int pairSign = t.charge * a.charge;
          if (cfgPairCharge.value == 1 && pairSign <= 0) {
            continue; // keep only like-sign pairs
          }
          if (cfgPairCharge.value == -1 && pairSign >= 0) {
            continue; // keep only unlike-sign pairs
          }
        }
        float deltaPhi = RecoDecay::constrainAngle(t.phi - a.phi, -o2::constants::math::PIHalf);
        float deltaEta = t.eta - a.eta;
        pairHist->Fill(step,
                       deltaEta,
                       a.pt,
                       t.pt,
                       multiplicity,
                       deltaPhi,
                       posZ,
                       eventWeight);
      }
    }
  }
  template <class CollType, class TCFTracks>
  void processSameDerivedPIDT(CollType const& collision, TCFTracks const& tracks, TracksPID const& tracksAll)
  {
    if (cfgVerbosity > 0) {
      LOGF(info, "processSameDerivedPIDT: Tracks for collision: %d | Vertex: %.1f | Multiplicity/Centrality: %.1f", tracks.size(), collision.posZ(), collision.multiplicity());
    }

    loadEfficiency(collision.timestamp());

    const auto multiplicity = eventMultiplicityFT0MOrFallback(collision, collision.multiplicity());

    using BinningTypeDerived = ColumnBinningPolicy<aod::collision::PosZ, aod::cfcollision::Multiplicity>;
    BinningTypeDerived configurableBinningDerived{{axisVertex, axisMultiplicity}, true};
    int bin = configurableBinningDerived.getBin({collision.posZ(), multiplicity});
    registry.fill(HIST("eventcount_same"), bin);
    registry.fill(HIST("trackcount_same"), bin, tracks.size());

    // Kinematic QA
    fillQA(collision, multiplicity, tracks);

    // PID-selected trigger/associate lists via CFTrackRefs -> AO2D TracksPID
    std::vector<SimpleTrack> triggerTracks;
    std::vector<SimpleTrack> associatedTracks;
    triggerTracks.reserve(tracks.size());
    associatedTracks.reserve(tracks.size());

    for (auto const& cftrk : tracks) {
      const auto trk = trackAt(tracksAll, cftrk.trackId());
      if (trk.sign() == 0) {
        continue;
      }

      if (passPIDForSpecies(trk, cfgTriggerSpecies.value)) {
        if (cfgTriggerCharge.value == 0 || trk.sign() == cfgTriggerCharge.value) {
          triggerTracks.push_back(SimpleTrack{cftrk.eta(), cftrk.phi(), cftrk.pt(), trk.sign()});
        }
      }

      if (passPIDForSpecies(trk, cfgAssociatedSpecies.value)) {
        if (cfgAssociatedCharge.value == 0 || trk.sign() == cfgAssociatedCharge.value) {
          associatedTracks.push_back(SimpleTrack{cftrk.eta(), cftrk.phi(), cftrk.pt(), trk.sign()});
        }
      }
    }

    if (triggerTracks.empty() || associatedTracks.empty()) {
      return;
    }

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsSimple(same,
                           CorrelationContainer::kCFStepReconstructed,
                           triggerTracks,
                           associatedTracks,
                           multiplicity,
                           collision.posZ(),
                           1.0f);

    if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
      same->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelationsSimple(same,
                             CorrelationContainer::kCFStepCorrected,
                             triggerTracks,
                             associatedTracks,
                             multiplicity,
                             collision.posZ(),
                             1.0f);
    }
  }

  void processSameDerivedPID(DerivedCollisions::iterator const& collision,
                             DerivedTracksWithRefs const& tracks,
                             TracksPID const& tracksAll)
  {
    processSameDerivedPIDT(collision, tracks, tracksAll);
  }
  PROCESS_SWITCH(Nucleibalance, processSameDerivedPID, "Process same event on derived data with PID via CFTrackRefs", false);

  void processSameDerivedMultSetPID(soa::Filtered<soa::Join<aod::CFCollisions, aod::CFMultSets>>::iterator const& collision,
                                    DerivedTracksWithRefs const& tracks,
                                    TracksPID const& tracksAll)
  {
    processSameDerivedPIDT(collision, tracks, tracksAll);
  }
  PROCESS_SWITCH(Nucleibalance, processSameDerivedMultSetPID, "Process same event on derived data with multiplicity sets and PID via CFTrackRefs", false);
  // AO2D-based processing: same + mixed events with PID-selected pi, p, d
  void processAO2D(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>::iterator const& collision,
                   TracksPID const& tracksAll)
  {
    // Event selection on vertex
    if (std::abs(collision.posZ()) > cfgCutVertex) {
      return;
    }

    // Event selection (cfgTrigger) -- AO2D only
    if (!keepCollisionAO2D(collision)) {
      return;
    }

    const auto thisCollIndex = collision.globalIndex();

    // Per-event containers
    std::vector<SimpleTrack> eventTracks;      // all selected tracks for multiplicity / QA
    std::vector<SimpleTrack> triggerTracks;    // PID + charge selected triggers
    std::vector<SimpleTrack> associatedTracks; // PID + charge selected associates

    // Loop over all tracks and select those belonging to this collision
    for (auto const& trk : tracksAll) {
      if (trk.collisionId() != thisCollIndex) {
        continue;
      }

      // Kinematic cuts
      if (std::abs(trk.eta()) > cfgCutEta || trk.pt() < cfgCutPt) {
        continue;
      }

      // Track-quality cuts (default: global tracks)
      if (!passTrackQualityAO2D(trk)) {
        continue;
      }

      // Save for multiplicity / QA (keep charge even if neutral)
      eventTracks.push_back(SimpleTrack{trk.eta(), trk.phi(), trk.pt(), static_cast<int>(trk.sign())});

      if (trk.sign() == 0) {
        continue;
      }

      // Trigger selection: PID + charge
      if (passPIDForSpecies(trk, cfgTriggerSpecies.value)) {
        if (cfgTriggerCharge.value == 0 || trk.sign() == cfgTriggerCharge.value) {
          triggerTracks.push_back(SimpleTrack{trk.eta(), trk.phi(), trk.pt(), trk.sign()});
        }
      }

      // Associated selection: PID + charge
      if (passPIDForSpecies(trk, cfgAssociatedSpecies.value)) {
        if (cfgAssociatedCharge.value == 0 || trk.sign() == cfgAssociatedCharge.value) {
          associatedTracks.push_back(SimpleTrack{trk.eta(), trk.phi(), trk.pt(), trk.sign()});
        }
      }
    }

    if (triggerTracks.empty() || associatedTracks.empty()) {
      return;
    }

    const float multiplicity =
      eventMultiplicityFT0MOrFallback(collision, static_cast<float>(eventTracks.size()));

    // QA on tracks for this event (AO2D-based)
    registry.fill(HIST("multiplicity"), multiplicity);
    for (const auto& t : eventTracks) {
      registry.fill(HIST("yields"), multiplicity, t.pt, t.eta);
      registry.fill(HIST("etaphi"), multiplicity, t.eta, t.phi);
    }

    // --------------------------
    // SAME-EVENT CORRELATIONS
    // --------------------------
    same->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsSimple(same,
                           CorrelationContainer::kCFStepReconstructed,
                           triggerTracks,
                           associatedTracks,
                           multiplicity,
                           collision.posZ(),
                           1.0f);

    // Optional efficiency-corrected step (if you configure efficiencies)
    if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
      if (!cfg.efficiencyLoaded) {
        // For AO2D only, the timestamp is not crucial if you use local efficiency files
        loadEfficiency(0);
      }
      same->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelationsSimple(same,
                             CorrelationContainer::kCFStepCorrected,
                             triggerTracks,
                             associatedTracks,
                             multiplicity,
                             collision.posZ(),
                             1.0f);
    }

    // --------------------------
    // MIXED-EVENT CORRELATIONS
    // --------------------------
    for (auto const& prev : mMixEvents) {
      if (prev.triggerTracks.empty() || prev.associatedTracks.empty()) {
        continue;
      }

      mixed->fillEvent(prev.multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsSimple(mixed,
                             CorrelationContainer::kCFStepReconstructed,
                             triggerTracks,
                             prev.associatedTracks,
                             prev.multiplicity,
                             prev.zvtx,
                             1.0f);

      if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
        mixed->fillEvent(prev.multiplicity, CorrelationContainer::kCFStepCorrected);
        fillCorrelationsSimple(mixed,
                               CorrelationContainer::kCFStepCorrected,
                               triggerTracks,
                               prev.associatedTracks,
                               prev.multiplicity,
                               prev.zvtx,
                               1.0f);
      }
    }

    // --------------------------
    // UPDATE MIXING BUFFER
    // --------------------------
    MixEventEntry entry;
    entry.multiplicity = multiplicity;
    entry.zvtx = collision.posZ();
    entry.triggerTracks = std::move(triggerTracks);
    entry.associatedTracks = std::move(associatedTracks);

    mMixEvents.push_front(std::move(entry));
    if (mMixEvents.size() > static_cast<size_t>(cfgNoMixedEvents.value)) {
      mMixEvents.pop_back();
    }
  }

  PROCESS_SWITCH(Nucleibalance, processAO2D,
                 "Process AO2D: same- and mixed-event correlations with PID", true);

  // AO2D-MC processing: reconstructed correlations + MC efficiency (truth vs reco)
  void processAO2DMC(CollisionsAO2DMC::iterator const& collision,
                     TracksPIDMC const& tracksAll,
                     aod::McCollisions const& /*mcCollisions*/,
                     aod::McParticles const& mcParticles)
  {
    if (std::abs(collision.posZ()) > cfgCutVertex) {
      return;
    }
    if (!keepCollisionAO2D(collision)) {
      return;
    }

    // Resolve MC collision index
    int mcCollIdx = -1;
    if constexpr (requires { collision.has_mcCollision(); }) {
      if (!collision.has_mcCollision()) {
        return;
      }
    }
    if constexpr (requires { collision.mcCollisionId(); }) {
      mcCollIdx = collision.mcCollisionId();
    } else if constexpr (requires { collision.mcCollision().globalIndex(); }) {
      mcCollIdx = collision.mcCollision().globalIndex();
    }
    if (mcCollIdx < 0) {
      return;
    }

    const auto thisCollIndex = collision.globalIndex();

    std::vector<SimpleTrack> eventTracks;
    std::vector<SimpleTrack> triggerTracks;
    std::vector<SimpleTrack> associatedTracks;

    for (auto const& trk : tracksAll) {
      if (trk.collisionId() != thisCollIndex) {
        continue;
      }
      if (std::abs(trk.eta()) > cfgCutEta || trk.pt() < cfgCutPt) {
        continue;
      }
      if (!passTrackQualityAO2D(trk)) {
        continue;
      }

      eventTracks.push_back(SimpleTrack{trk.eta(), trk.phi(), trk.pt(), static_cast<int>(trk.sign())});

      if (trk.sign() == 0) {
        continue;
      }

      if (passPIDForSpecies(trk, cfgTriggerSpecies.value) &&
          (cfgTriggerCharge.value == 0 || trk.sign() == cfgTriggerCharge.value)) {
        triggerTracks.push_back(SimpleTrack{trk.eta(), trk.phi(), trk.pt(), trk.sign()});
      }

      if (passPIDForSpecies(trk, cfgAssociatedSpecies.value) &&
          (cfgAssociatedCharge.value == 0 || trk.sign() == cfgAssociatedCharge.value)) {
        associatedTracks.push_back(SimpleTrack{trk.eta(), trk.phi(), trk.pt(), trk.sign()});
      }
    }

    if (triggerTracks.empty() || associatedTracks.empty()) {
      return;
    }

    const float multiplicity =
      eventMultiplicityFT0MOrFallback(collision, static_cast<float>(eventTracks.size()));

    // ---- MC efficiency filling (AO2D-MC) ----
    auto groupedMcParticles = mcParticles.sliceBy(mcParticlesPerCollision, mcCollIdx);

    for (auto const& mcPart : groupedMcParticles) {
      if (std::abs(mcPart.eta()) > cfgCutMCEta || mcPart.pt() < cfgCutMCPt) {
        continue;
      }
      if (mcPart.isPhysicalPrimary() && chargeFromPdg(mcPart.pdgCode()) != 0) {
        same->getTrackHistEfficiency()->Fill(CorrelationContainer::MC,
                                             mcPart.eta(), mcPart.pt(),
                                             getSpecies(mcPart.pdgCode()),
                                             multiplicity, collision.posZ());
      }
    }

    for (auto const& trk : tracksAll) {
      if (trk.collisionId() != thisCollIndex) {
        continue;
      }
      if (std::abs(trk.eta()) > cfgCutEta || trk.pt() < cfgCutPt) {
        continue;
      }
      if (!passTrackQualityAO2D(trk)) {
        continue;
      }

      if constexpr (requires { trk.has_mcParticle(); trk.mcParticle(); }) {
        if (trk.has_mcParticle()) {
          const auto mcPart = trk.mcParticle();
          if (mcPart.isPhysicalPrimary()) {
            same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoPrimaries,
                                                 mcPart.eta(), mcPart.pt(),
                                                 getSpecies(mcPart.pdgCode()),
                                                 multiplicity, collision.posZ());
          }
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoAll,
                                               mcPart.eta(), mcPart.pt(),
                                               getSpecies(mcPart.pdgCode()),
                                               multiplicity, collision.posZ());
        } else {
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::Fake,
                                               trk.eta(), trk.pt(),
                                               0, multiplicity, collision.posZ());
        }
      }
    }

    // ---- QA ----
    registry.fill(HIST("multiplicity"), multiplicity);
    for (const auto& t : eventTracks) {
      registry.fill(HIST("yields"), multiplicity, t.pt, t.eta);
      registry.fill(HIST("etaphi"), multiplicity, t.eta, t.phi);
    }

    // ---- Same-event ----
    same->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsSimple(same, CorrelationContainer::kCFStepReconstructed,
                           triggerTracks, associatedTracks,
                           multiplicity, collision.posZ(), 1.0f);

    if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
      if (!cfg.efficiencyLoaded) {
        loadEfficiency(0);
      }
      same->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelationsSimple(same, CorrelationContainer::kCFStepCorrected,
                             triggerTracks, associatedTracks,
                             multiplicity, collision.posZ(), 1.0f);
    }

    // ---- Mixed-event ----
    for (auto const& prev : mMixEventsMC) {
      if (prev.triggerTracks.empty() || prev.associatedTracks.empty()) {
        continue;
      }

      mixed->fillEvent(prev.multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsSimple(mixed, CorrelationContainer::kCFStepReconstructed,
                             triggerTracks, prev.associatedTracks,
                             prev.multiplicity, prev.zvtx, 1.0f);

      if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
        mixed->fillEvent(prev.multiplicity, CorrelationContainer::kCFStepCorrected);
        fillCorrelationsSimple(mixed, CorrelationContainer::kCFStepCorrected,
                               triggerTracks, prev.associatedTracks,
                               prev.multiplicity, prev.zvtx, 1.0f);
      }
    }

    MixEventEntry entry;
    entry.multiplicity = multiplicity;
    entry.zvtx = collision.posZ();
    entry.triggerTracks = std::move(triggerTracks);
    entry.associatedTracks = std::move(associatedTracks);

    mMixEventsMC.push_front(std::move(entry));
    if (mMixEventsMC.size() > static_cast<size_t>(cfgNoMixedEvents.value)) {
      mMixEventsMC.pop_back();
    }
  }

  PROCESS_SWITCH(Nucleibalance, processAO2DMC,
                 "Process AO2D-MC: reconstructed correlations + MC efficiency via labels", false);

  template <class CollType, class TTracks1, class TTracks2>
  void processSameDerivedT(CollType const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    using BinningTypeDerived = ColumnBinningPolicy<aod::collision::PosZ, aod::cfcollision::Multiplicity>;
    BinningTypeDerived configurableBinningDerived{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    if (cfgVerbosity > 0) {
      LOGF(info, "processSameDerivedT: Tracks for collision: %d/%d | Vertex: %.1f | Multiplicity/Centrality: %.1f", tracks1.size(), tracks2.size(), collision.posZ(), collision.multiplicity());
    }
    loadEfficiency(collision.timestamp());

    const auto multiplicity = eventMultiplicityFT0MOrFallback(collision, collision.multiplicity());
    int field = 0;
    if (cfgTwoTrackCut > 0) {
      field = getMagneticField(collision.timestamp());
    }

    int bin = configurableBinningDerived.getBin({collision.posZ(), multiplicity});
    registry.fill(HIST("eventcount_same"), bin);
    registry.fill(HIST("trackcount_same"), bin, tracks1.size());

    fillQA(collision, multiplicity, tracks1);

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f);

    if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
      same->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelations<CorrelationContainer::kCFStepCorrected>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f);
    }
  }

  void processSameDerived(DerivedCollisions::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    processSameDerivedT(collision, tracks, tracks);
  }
  PROCESS_SWITCH(Nucleibalance, processSameDerived, "Process same event on derived data", false);

  void processSameDerivedMultSet(soa::Filtered<soa::Join<aod::CFCollisions, aod::CFMultSets>>::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    processSameDerivedT(collision, tracks, tracks);
  }
  PROCESS_SWITCH(Nucleibalance, processSameDerivedMultSet, "Process same event on derived data with multiplicity sets", false);

  template <class CollType, typename... TrackTypes>
  void processMixedDerivedT(CollType const& collisions, TrackTypes&&... tracks)
  {
    auto getMultiplicity =
      [this](auto& col) {
        (void)this; // fix compile error on unused 'this' capture
        return eventMultiplicityFT0MOrFallback(col, col.multiplicity());
      };

    using BinningTypeDerived = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    BinningTypeDerived configurableBinningDerived{{getMultiplicity}, {axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    //  Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TA = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TB = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<CollType, TA, TB, BinningTypeDerived> pairs{configurableBinningDerived, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      float multiplicity = getMultiplicity(collision1);
      int bin = configurableBinningDerived.getBin(std::tuple(collision1.posZ(), multiplicity));
      float eventWeight = 1.0f / it.currentWindowNeighbours();
      int field = 0;
      if (cfgTwoTrackCut > 0) {
        field = getMagneticField(collision1.timestamp());
      }

      if (cfgVerbosity > 0) {
        LOGF(info, "processMixedDerived: Mixed collisions bin: %d pair: [%d, %d] %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, it.isNewWindow(), it.currentWindowNeighbours(), collision1.globalIndex(), collision1.posZ(), collision1.multiplicity(), collision2.globalIndex(), collision2.posZ(), collision2.multiplicity());
      }

      if (it.isNewWindow()) {
        loadEfficiency(collision1.timestamp());
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      }

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      registry.fill(HIST("eventcount_mixed"), bin);
      registry.fill(HIST("trackcount_mixed"), bin, tracks1.size(), tracks2.size());
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), field, eventWeight);

      if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
        if (it.isNewWindow()) {
          mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
        }
        fillCorrelations<CorrelationContainer::kCFStepCorrected>(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), field, eventWeight);
      }
    }
  }

  void processMixedDerived(DerivedCollisions const& collisions, DerivedTracks const& tracks)
  {
    processMixedDerivedT(collisions, tracks);
  }
  PROCESS_SWITCH(Nucleibalance, processMixedDerived, "Process mixed events on derived data", false);

  void processMixedDerivedMultSet(soa::Filtered<soa::Join<aod::CFCollisions, aod::CFMultSets>> const& collisions, DerivedTracks const& tracks)
  {
    processMixedDerivedT(collisions, tracks);
  }
  PROCESS_SWITCH(Nucleibalance, processMixedDerivedMultSet, "Process mixed events on derived data with multiplicity sets", false);

  // Mixed-event processing on derived data with PID via CFTrackRefs -> AO2D TracksPID
  template <class CollType, typename... TrackTypes>
  void processMixedDerivedPIDT(CollType const& collisions, TracksPID const& tracksAll, TrackTypes&&... tracks)
  {
    auto getMultiplicity =
      [this](auto& col) {
        (void)this;
        return eventMultiplicityFT0MOrFallback(col, col.multiplicity());
      };

    using BinningTypeDerived = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    BinningTypeDerived configurableBinningDerived{{getMultiplicity}, {axisVertex, axisMultiplicity}, true};

    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TA = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TB = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<CollType, TA, TB, BinningTypeDerived> pairs{configurableBinningDerived, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};

    for (auto it = pairs.begin(); it != pairs.end(); ++it) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      const float multiplicity = getMultiplicity(collision1);
      const int bin = configurableBinningDerived.getBin(std::tuple(collision1.posZ(), multiplicity));
      const float eventWeight = 1.0f / it.currentWindowNeighbours();

      if (cfgVerbosity > 0) {
        LOGF(info, "processMixedDerivedPID: Mixed collisions bin: %d pair: [%d, %d] %d (%.3f, %.3f), %d (%.3f, %.3f)",
             bin, it.isNewWindow(), it.currentWindowNeighbours(),
             collision1.globalIndex(), collision1.posZ(), collision1.multiplicity(),
             collision2.globalIndex(), collision2.posZ(), collision2.multiplicity());
      }

      if (it.isNewWindow()) {
        loadEfficiency(collision1.timestamp());
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      }

      registry.fill(HIST("eventcount_mixed"), bin);
      registry.fill(HIST("trackcount_mixed"), bin, tracks1.size(), tracks2.size());

      std::vector<SimpleTrack> triggerTracks;
      std::vector<SimpleTrack> associatedTracks;
      triggerTracks.reserve(tracks1.size());
      associatedTracks.reserve(tracks2.size());

      // Triggers from collision1
      for (auto const& cftrk : tracks1) {
        const auto trk = trackAt(tracksAll, cftrk.trackId());
        if (trk.sign() == 0) {
          continue;
        }
        if (passPIDForSpecies(trk, cfgTriggerSpecies.value)) {
          if (cfgTriggerCharge.value == 0 || trk.sign() == cfgTriggerCharge.value) {
            triggerTracks.push_back(SimpleTrack{cftrk.eta(), cftrk.phi(), cftrk.pt(), trk.sign()});
          }
        }
      }

      // Associates from collision2
      for (auto const& cftrk : tracks2) {
        const auto trk = trackAt(tracksAll, cftrk.trackId());
        if (trk.sign() == 0) {
          continue;
        }
        if (passPIDForSpecies(trk, cfgAssociatedSpecies.value)) {
          if (cfgAssociatedCharge.value == 0 || trk.sign() == cfgAssociatedCharge.value) {
            associatedTracks.push_back(SimpleTrack{cftrk.eta(), cftrk.phi(), cftrk.pt(), trk.sign()});
          }
        }
      }

      if (!triggerTracks.empty() && !associatedTracks.empty()) {
        fillCorrelationsSimple(mixed,
                               CorrelationContainer::kCFStepReconstructed,
                               triggerTracks,
                               associatedTracks,
                               multiplicity,
                               collision1.posZ(),
                               eventWeight);
      }

      if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
        if (it.isNewWindow()) {
          mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
        }
        if (!triggerTracks.empty() && !associatedTracks.empty()) {
          fillCorrelationsSimple(mixed,
                                 CorrelationContainer::kCFStepCorrected,
                                 triggerTracks,
                                 associatedTracks,
                                 multiplicity,
                                 collision1.posZ(),
                                 eventWeight);
        }
      }
    }
  }

  void processMixedDerivedPID(DerivedCollisions const& collisions,
                              DerivedTracksWithRefs const& tracks,
                              TracksPID const& tracksAll)
  {
    processMixedDerivedPIDT(collisions, tracksAll, tracks);
  }
  PROCESS_SWITCH(Nucleibalance, processMixedDerivedPID, "Process mixed events on derived data with PID via CFTrackRefs", false);

  void processMixedDerivedMultSetPID(soa::Filtered<soa::Join<aod::CFCollisions, aod::CFMultSets>> const& collisions,
                                     DerivedTracksWithRefs const& tracks,
                                     TracksPID const& tracksAll)
  {
    processMixedDerivedPIDT(collisions, tracksAll, tracks);
  }
  PROCESS_SWITCH(Nucleibalance, processMixedDerivedMultSetPID, "Process mixed events on derived data with multiplicity sets and PID via CFTrackRefs", false);

  int getSpecies(int pdgCode)
  {
    switch (pdgCode) {
      case PdgPion: // pion
      case -PdgPion:
        return 0;
      case PdgKaon: // Kaon
      case -PdgKaon:
        return 1;
      case PdgProton: // proton
      case -PdgProton:
        return 2;
    }
    if (std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), pdgCode) != cfgMcTriggerPDGs->end())
      return 4;
    else
      return 3;
  }

  // NOTE SmallGroups includes soa::Filtered always
  Preslice<aod::CFTracksWithLabel> perCollision = aod::cftrack::cfCollisionId;
  void processMCEfficiency(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, aod::CFMcParticles const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions, aod::CFTracksWithLabel const& tracks)
  {
    if (cfgVerbosity > 0) {
      LOGF(info, "MC collision at vtx-z = %f with %d mc particles and %d reconstructed collisions", mcCollision.posZ(), mcParticles.size(), collisions.size());
    }

    auto multiplicity = mcCollision.multiplicity();
    if (cfgCentBinsForMC > 0) {
      if (collisions.size() == 0) {
        return;
      }
      for (const auto& collision : collisions) {
        multiplicity = collision.multiplicity();
      }
    }
    // Primaries
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() && mcParticle.sign() != 0) {
        same->getTrackHistEfficiency()->Fill(CorrelationContainer::MC, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
      }
    }
    for (const auto& collision : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "  Reconstructed collision at vtx-z = %f", collision.posZ());
        LOGF(info, "  which has %d tracks", groupedTracks.size());
      }

      for (const auto& track : groupedTracks) {
        if (track.has_cfMCParticle()) {
          const auto& mcParticle = track.cfMCParticle();
          if (mcParticle.isPhysicalPrimary()) {
            same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoPrimaries, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
          }
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoAll, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
          // LOGF(info, "Filled track %d", track.globalIndex());
        } else {
          // fake track
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::Fake, track.eta(), track.pt(), 0, multiplicity, mcCollision.posZ());
        }
      }
    }
  }
  PROCESS_SWITCH(Nucleibalance, processMCEfficiency, "MC: Extract efficiencies", false);

  template <class Particles1, class Particles2>
  void processMCSameDerivedT(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, Particles1 const& mcParticles1, Particles2 const& mcParticles2, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions)
  {
    if (cfgVerbosity > 0) {
      LOGF(info, "processMCSameDerivedT. MC collision: %d, particles1: %d, particles2: %d, collisions: %d", mcCollision.globalIndex(), mcParticles1.size(), mcParticles2.size(), collisions.size());
    }

    auto multiplicity = mcCollision.multiplicity();
    if (cfgCentBinsForMC > 0) {
      if (collisions.size() == 0) {
        return;
      }
      for (const auto& collision : collisions) {
        multiplicity = collision.multiplicity();
      }
    }

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    fillCorrelations<CorrelationContainer::kCFStepAll>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    if (collisions.size() == 0) {
      return;
    }

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepVertex);
    fillCorrelations<CorrelationContainer::kCFStepVertex>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepTrackedOnlyPrim);
    fillCorrelations<CorrelationContainer::kCFStepTrackedOnlyPrim>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepTracked);
    fillCorrelations<CorrelationContainer::kCFStepTracked>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    // NOTE kCFStepReconstructed and kCFStepCorrected are filled in processSameDerived
    //      This also means that if a MC collision had several reconstructed vertices (collisions), all of them are filled
  }

  // NOTE SmallGroups includes soa::Filtered always
  void processMCSameDerived(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, soa::Filtered<aod::CFMcParticles> const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions) // TODO. For mixed no need to check the daughters since the events are different
  {
    processMCSameDerivedT(mcCollision, mcParticles, mcParticles, collisions);
  }
  PROCESS_SWITCH(Nucleibalance, processMCSameDerived, "Process MC same event on derived data", false);

  PresliceUnsorted<aod::CFCollisionsWithLabel> collisionPerMCCollision = aod::cfcollision::cfMcCollisionId;
  template <typename... ParticleTypes>
  void processMCMixedDerivedT(soa::Filtered<aod::CFMcCollisions> const& mcCollisions, soa::Filtered<aod::CFCollisionsWithLabel> const& collisions, ParticleTypes&&... particles)
  {
    bool useMCMultiplicity = (cfgCentBinsForMC == 0);
    auto getMultiplicity =
      [&collisions, &useMCMultiplicity, this](auto& col) {
        if (useMCMultiplicity)
          return col.multiplicity();
        auto groupedCollisions = collisions.sliceBy(collisionPerMCCollision, col.globalIndex());
        if (groupedCollisions.size() == 0)
          return -1.0f;
        return groupedCollisions.begin().multiplicity();
      };

    using BinningTypeMCDerived = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::mccollision::PosZ, decltype(getMultiplicity)>;
    BinningTypeMCDerived configurableBinning{{getMultiplicity}, {axisVertex, axisMultiplicity}, true};

    // Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    auto tuple = std::make_tuple(std::forward<ParticleTypes>(particles)...);
    using TA = std::tuple_element<0, decltype(tuple)>::type;
    using TB = std::tuple_element<std::tuple_size_v<decltype(tuple)> - 1, decltype(tuple)>::type;
    Pair<soa::Filtered<aod::CFMcCollisions>, TA, TB, BinningTypeMCDerived> pairs{configurableBinning, cfgNoMixedEvents, -1, mcCollisions, tuple, &cache}; // -1 is the number of the bin to skip

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      float eventWeight = 1.0f / it.currentWindowNeighbours();

      float multiplicity = getMultiplicity(collision1);
      if (cfgVerbosity > 0) {
        int bin = configurableBinning.getBin(std::tuple(collision1.posZ(), multiplicity));
        LOGF(info, "processMCMixedDerived: Mixed collisions bin: %d pair: [%d, %d] %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, it.isNewWindow(), it.currentWindowNeighbours(), collision1.globalIndex(), collision1.posZ(), getMultiplicity(collision1), collision2.globalIndex(), collision2.posZ(), getMultiplicity(collision2));
      }

      // STEP 0
      if (it.isNewWindow()) {
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
      }
      fillCorrelations<CorrelationContainer::kCFStepAll>(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), 0, eventWeight);
      // check if collision1 has at least one reconstructed collision
      auto groupedCollisions = collisions.sliceBy(collisionPerMCCollision, collision1.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "Found %d related collisions", groupedCollisions.size());
      }
      if (groupedCollisions.size() == 0) {
        continue;
      }

      // STEP 2, 4, 5
      if (it.isNewWindow()) {
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepVertex);
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepTrackedOnlyPrim);
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepTracked);
      }
      fillCorrelations<CorrelationContainer::kCFStepVertex>(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), 0, eventWeight);
      fillCorrelations<CorrelationContainer::kCFStepTrackedOnlyPrim>(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), 0, eventWeight);
      fillCorrelations<CorrelationContainer::kCFStepTracked>(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), 0, eventWeight);

      // NOTE kCFStepReconstructed and kCFStepCorrected are filled in processMixedDerived
      //      This also means that if a MC collision had several reconstructed vertices (collisions), all of them are filled
    }
  }

  void processMCMixedDerived(soa::Filtered<aod::CFMcCollisions> const& mcCollisions, soa::Filtered<aod::CFMcParticles> const& mcParticles, soa::Filtered<aod::CFCollisionsWithLabel> const& collisions)
  {
    processMCMixedDerivedT(mcCollisions, collisions, mcParticles);
  }
  PROCESS_SWITCH(Nucleibalance, processMCMixedDerived, "Process MC mixed events on derived data", false);
};

// Lambda* proxy analysis task based on deuteron proxy
struct Lambdastarproxy {
  // ---- Ion/nucleus PDG encoding helpers (10LZZZAAAI) ----
  // Note: these are *format* constants (not particle PDG species codes)
  static constexpr int IonCodeThreshold = 1000000000; // 10^9
  static constexpr int IonZDivisor = 10000;
  static constexpr int IonZModulo = 1000;
  // ---- Common PDG codes (use named values; avoid explicit literals) ----
  static constexpr int PdgElectron = static_cast<int>(PDG_t::kElectron);
  static constexpr int PdgMuon = static_cast<int>(PDG_t::kMuonMinus);
  static constexpr int PdgPion = static_cast<int>(PDG_t::kPiPlus);
  static constexpr int PdgKaon = static_cast<int>(PDG_t::kKPlus);
  static constexpr int PdgProton = static_cast<int>(PDG_t::kProton);
  // ---- Named defaults to avoid magic numbers (o2-linter) ----
  static constexpr float CutVertexDefault = 10.f;
  static constexpr float CutPtMinDefault = 0.5f;
  static constexpr float CutEtaMaxDefault = 0.8f;

  static constexpr float CutMCPtMinDefault = 0.1f;
  static constexpr float CutMCEtaMaxDefault = 0.8f;
  static constexpr int FillMCTruthDefault = 1;

  static constexpr float NsigmaTPCDefault = 3.f;
  static constexpr float NsigmaTOFDefault = 3.f;

  static constexpr bool RequireGlobalTrackDefault = true;
  static constexpr int TPCNClsMinDefault = 70;
  static constexpr float DcaXYMaxDefault = 0.1f;
  static constexpr float DcaZMaxDefault = 0.2f;
  static constexpr float Chi2PerTPCClusterDefault = 2.5f;
  static constexpr float Chi2PerITSClusterDefault = 36.f;
  static constexpr int ITSNClustersDefault = 5;

  static constexpr int TriggerDefault = 4;
  // Named trigger codes
  static constexpr int TriggerNone = 0;
  static constexpr int TriggerSel8 = 1;
  static constexpr int TriggerSel8Quality = 2;
  static constexpr int TriggerSel8OccQuality = 3;
  static constexpr int TriggerSel8NoSbpZvtx = 4;
  static constexpr int MinOccDefault = 0;
  static constexpr int MaxOccDefault = 3000;

  static constexpr int NoMixedEventsDefault = 5;
  static constexpr float MixZvtxMaxDefault = 2.0f;
  static constexpr float MixMultMaxDefault = 50.0f;

  static constexpr float ProxyMomentumScale = 0.5f;
  static constexpr float TofBetaMin = 0.01f;
  static constexpr float TofBetaMax = 1.2f;
  static constexpr double Half = 0.5;
  // Basic configuration for event and track selection
  Configurable<float> lstarCutVertex{"lstarCutVertex", float{CutVertexDefault}, "Accepted z-vertex range (cm)"};
  Configurable<float> lstarCutPtMin{"lstarCutPtMin", float{CutPtMinDefault}, "Minimal pT for tracks (GeV/c)"};
  Configurable<float> lstarCutEtaMax{"lstarCutEtaMax", float{CutEtaMaxDefault}, "Max |eta| for tracks"};

  // MC truth
  Configurable<float> lstarCutMCPtMin{"lstarCutMCPtMin", float{CutMCPtMinDefault}, "Minimal pT for MC particles (GeV/c)"};
  Configurable<float> lstarCutMCEtaMax{"lstarCutMCEtaMax", float{CutMCEtaMaxDefault}, "Max |eta| for MC particles"};
  Configurable<int> lstarFillMCTruth{"lstarFillMCTruth", int{FillMCTruthDefault}, "Fill MC truth and reco-matching QA (AO2D-MC mode)"};

  // PID cuts
  Configurable<float> lstarCutNsigmaTPCPi{"lstarCutNsigmaTPCPi", float{NsigmaTPCDefault}, "|nSigma^{TPC}_{#pi}| cut"};
  Configurable<float> lstarCutNsigmaTOFPi{"lstarCutNsigmaTOFPi", float{NsigmaTOFDefault}, "|nSigma^{TOF}_{#pi}| cut"};
  Configurable<float> lstarCutNsigmaTPCPr{"lstarCutNsigmaTPCPr", float{NsigmaTPCDefault}, "|nSigma^{TPC}_{p}| cut"};
  Configurable<float> lstarCutNsigmaTOFPr{"lstarCutNsigmaTOFPr", float{NsigmaTOFDefault}, "|nSigma^{TOF}_{p}| cut"};
  Configurable<float> lstarCutNsigmaTPCKaon{"lstarCutNsigmaTPCKaon", float{NsigmaTPCDefault}, "|nSigma^{TPC}_{K}| cut"};
  Configurable<float> lstarCutNsigmaTOFKaon{"lstarCutNsigmaTOFKaon", float{NsigmaTOFDefault}, "|nSigma^{TOF}_{K}| cut"};
  Configurable<float> lstarCutNsigmaTPCDe{"lstarCutNsigmaTPCDe", float{NsigmaTPCDefault}, "|nSigma^{TPC}_{d}| cut"};
  Configurable<float> lstarCutNsigmaTOFDe{"lstarCutNsigmaTOFDe", float{NsigmaTOFDefault}, "|nSigma^{TOF}_{d}| cut"};

  // Track quality
  Configurable<bool> lstarRequireGlobalTrack{"lstarRequireGlobalTrack", bool{RequireGlobalTrackDefault}, "Require global tracks (default)"};
  Configurable<int> lstarTPCNClsMin{"lstarTPCNClsMin", int{TPCNClsMinDefault}, "Minimum number of TPC clusters (tpcNClsFound)"};
  Configurable<float> lstarDcaXYMax{"lstarDcaXYMax", float{DcaXYMaxDefault}, "Max |DCA_{xy}| to PV (cm)"};
  Configurable<float> lstarDcaZMax{"lstarDcaZMax", float{DcaZMaxDefault}, "Max |DCA_{z}| to PV (cm)"};
  Configurable<float> lstarChi2PerTPCCluster{"lstarChi2PerTPCCluster", float{Chi2PerTPCClusterDefault}, "Maximum Chi2/cluster for the TPC track segment"};
  Configurable<float> lstarChi2PerITSCluster{"lstarChi2PerITSCluster", float{Chi2PerITSClusterDefault}, "Maximum Chi2/cluster for the ITS track segment"};
  Configurable<int> lstarITSNClusters{"lstarITSNClusters", int{ITSNClustersDefault}, "Minimum number of ITS clusters"};

  // Trigger + occupancy
  Configurable<int> lstarCfgTrigger{"lstarCfgTrigger", int{TriggerDefault}, "Event trigger selection: 0=none, 1=sel8, 2=sel8+(kNoSameBunchPileup+kIsGoodZvtxFT0vsPV+kIsGoodITSLayersAll), 3=sel8+occupancy+(kNoCollInTimeRangeStandard)+(kNoSameBunchPileup+kIsGoodZvtxFT0vsPV+kIsGoodITSLayersAll), 4=sel8+(kNoSameBunchPileup+kIsGoodZvtxFT0vsPV)"};
  Configurable<int> lstarMinOcc{"lstarMinOcc", int{MinOccDefault}, "minimum occupancy selection (for cfgTrigger==3)"};
  Configurable<int> lstarMaxOcc{"lstarMaxOcc", int{MaxOccDefault}, "maximum occupancy selection (for cfgTrigger==3)"};
  Configurable<int> lstarUseFT0M{"lstarUseFT0M", 1, "Use FT0M centrality (0-100) as multiplicity axis (1=ON, 0=use Ntracks)"};

  // --- Mixed-event configuration for pK / proxy invariant-mass background (AO2D only) ---
  Configurable<int> lstarNoMixedEvents{"lstarNoMixedEvents", int{NoMixedEventsDefault}, "Number of previous events kept for mixed-event background"};
  Configurable<float> lstarMixZvtxMax{"lstarMixZvtxMax", float{MixZvtxMaxDefault}, "Max |Δzvtx| (cm) for event mixing"};
  Configurable<float> lstarMixMultMax{"lstarMixMultMax", float{MixMultMaxDefault}, "Max |Δmult| for event mixing"};
  Configurable<int> lstarEnablePidQA{"lstarEnablePidQA", 0, "Enable PID QA histograms (dE/dx, TOF #beta, proxy invariant-mass QA, etc.): 1 = ON, 0 = OFF"};
  Configurable<int> lstarEnableSparse{"lstarEnableSparse", 0, "Enable THnSparse invariant-mass histograms (#Lambda^{*} pK and proxy); 1 = ON, 0 = OFF"};

  struct KaonCand {
    float px, py, pz;
    int charge;
    int tid;
  };
  struct ProxyCand {
    float px, py, pz;
    int charge;
    int tid;
  };

  // Helpers for invariant-mass kinematics
  static float phiFromPxPy(float px, float py)
  {
    return std::atan2(py, px);
  }

  static float rapidityFromEPz(double e, double pz)
  {
    const double num = e + pz;
    const double den = e - pz;
    if (num <= 0.0 || den <= 0.0) {
      return 0.f;
    }
    return static_cast<float>(Half * std::log(num / den));
  }

  // Mixed-event pool entry for pK / proxy background (AO2D only)
  struct LStarMixEventEntry {
    float mult = 0.f;
    float zvtx = 0.f;
    std::vector<KaonCand> kaons;
    std::vector<ProxyCand> proxies;
  };

  // Keep last N events for event-mixing
  std::deque<LStarMixEventEntry> mLStarMixEvents;

  template <typename TCollision>
  bool keepCollisionAO2D(TCollision const& collision) const
  {
    if (lstarCfgTrigger.value == TriggerNone) {
      return true;
    } else if (lstarCfgTrigger.value == TriggerSel8) {
      return collision.sel8();
    } else if (lstarCfgTrigger.value == TriggerSel8Quality) {
      return collision.sel8() &&
             collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
             collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
             collision.selection_bit(aod::evsel::kIsGoodITSLayersAll);
    } else if (lstarCfgTrigger.value == TriggerSel8OccQuality) {
      const int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < lstarMinOcc.value || occupancy >= lstarMaxOcc.value) {
        return false;
      }
      return collision.sel8() &&
             collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
             collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
             collision.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) &&
             collision.selection_bit(aod::evsel::kIsGoodITSLayersAll);
    } else if (lstarCfgTrigger.value == TriggerSel8NoSbpZvtx) {
      return collision.sel8() &&
             collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
             collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV);
    }

    LOGF(warn, "Invalid lstarCfgTrigger=%d. Accepting all collisions.", lstarCfgTrigger.value);
    return true;
  }

  template <typename TCollision>
  float eventMultiplicityFT0MOrFallback(TCollision const& collision, float fallback) const
  {
    if (lstarUseFT0M.value == 0) {
      return fallback;
    }

    if constexpr (requires { collision.centFT0M(); }) {
      const float v = collision.centFT0M();
      if (v >= 0.f) {
        return v;
      }
    }

    if constexpr (requires { collision.centFT0MValid(); }) {
      if (collision.centFT0MValid()) {
        if constexpr (requires { collision.centFT0M(); }) {
          return collision.centFT0M();
        }
      }
    }

    return fallback;
  }

  static int chargeFromPdg(int pdg)
  {
    const int apdg = std::abs(pdg);

    // Ions/nuclei: PDG code format 10LZZZAAAI -> Z is encoded in digits [7..5]
    if (apdg >= IonCodeThreshold) {
      const int z = (apdg / IonZDivisor) % IonZModulo;
      return (pdg >= 0) ? z : -z;
    }

    switch (apdg) {
      case PdgElectron: // e
      case PdgMuon:     // mu
      case PdgPion:     // pi
      case PdgKaon:     // K
      case PdgProton:   // p
        return (pdg >= 0) ? 1 : -1;
      default:
        return 0;
    }
  }

  // Histogram registry for this task
  HistogramRegistry histos{"lstarRegistry"};

  // Filters
  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < lstarCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < lstarCutEtaMax) && (aod::track::pt > lstarCutPtMin);

  // Tracks with PID information from TPC and TOF services for all relevant species
  // NOTE: aod::TracksExtra is needed for TPC dE/dx (tpcSignal)
  using TracksWithPID = soa::Join<aod::Tracks,
                                  aod::TracksExtra,
                                  aod::TrackSelection,
                                  aod::TracksDCA,
                                  aod::pidTPCFullPi,
                                  aod::pidTOFFullPi,
                                  aod::pidTPCFullKa,
                                  aod::pidTOFFullKa,
                                  aod::pidTPCFullPr,
                                  aod::pidTOFFullPr,
                                  aod::pidTPCFullDe,
                                  aod::pidTOFFullDe,
                                  aod::pidTOFbeta>;
  template <typename TTrack>
  bool passTrackQuality(const TTrack& trk) const
  {
    if (lstarRequireGlobalTrack.value) {
      if constexpr (requires { trk.isGlobalTrack(); }) {
        if (!trk.isGlobalTrack()) {
          return false;
        }
      } else if constexpr (requires { trk.isGlobalTrackSDD(); }) {
        if (!trk.isGlobalTrackSDD()) {
          return false;
        }
      }
    }

    if constexpr (requires { trk.itsNCls(); }) {
      if (lstarITSNClusters.value > 0 && trk.itsNCls() < lstarITSNClusters.value) {
        return false;
      }
    }

    if constexpr (requires { trk.tpcNClsFound(); }) {
      if (lstarTPCNClsMin.value > 0 && trk.tpcNClsFound() < lstarTPCNClsMin.value) {
        return false;
      }
    }

    if constexpr (requires { trk.tpcChi2NCl(); }) {
      if (lstarChi2PerTPCCluster.value > 0.f && trk.tpcChi2NCl() > lstarChi2PerTPCCluster.value) {
        return false;
      }
    }

    if constexpr (requires { trk.itsChi2NCl(); }) {
      if (lstarChi2PerITSCluster.value > 0.f && trk.itsChi2NCl() > lstarChi2PerITSCluster.value) {
        return false;
      }
    }

    if constexpr (requires { trk.dcaXY(); }) {
      if (lstarDcaXYMax.value > 0.f && std::abs(trk.dcaXY()) > lstarDcaXYMax.value) {
        return false;
      }
    }

    if constexpr (requires { trk.dcaZ(); }) {
      if (lstarDcaZMax.value > 0.f && std::abs(trk.dcaZ()) > lstarDcaZMax.value) {
        return false;
      }
    }

    return true;
  }
  using CollisionsWithEvSel = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using FilteredCollisions = soa::Filtered<CollisionsWithEvSel>;
  using FilteredTracks = soa::Filtered<TracksWithPID>;

  // AO2D-MC variants (for truth QA and reco->MC matching)
  using CollisionsWithEvSelMC = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::McCollisionLabels>;
  using FilteredCollisionsMC = soa::Filtered<CollisionsWithEvSelMC>;

  using TracksWithPIDMC = soa::Join<aod::Tracks,
                                    aod::TracksExtra,
                                    aod::TrackSelection,
                                    aod::TracksDCA,
                                    aod::pidTPCFullPi,
                                    aod::pidTOFFullPi,
                                    aod::pidTPCFullKa,
                                    aod::pidTOFFullKa,
                                    aod::pidTPCFullPr,
                                    aod::pidTOFFullPr,
                                    aod::pidTPCFullDe,
                                    aod::pidTOFFullDe,
                                    aod::pidTOFbeta,
                                    aod::McTrackLabels>;
  using FilteredTracksMC = soa::Filtered<TracksWithPIDMC>;

  Preslice<aod::McParticles> mcParticlesPerCollision = aod::mcparticle::mcCollisionId;

  void init(o2::framework::InitContext&)
  {
    AxisSpec massAxis{200, 1.4, 1.9, "M_{pK} (GeV/c^{2})"};
    AxisSpec ptAxis{100, 0., 10., "p_{T} (GeV/c)"};
    AxisSpec nsAxis{100, -10., 10., "n#sigma"};
    AxisSpec pAxis{100, 0., 10., "p (GeV/c)"};
    AxisSpec etaAxis{80, -2., 2., "#eta"};
    AxisSpec phiAxis{64, 0., o2::constants::math::TwoPI, "#varphi"};
    AxisSpec centAxis{100, 0., 100., "multiplicity/centrality"};

    AxisSpec pdgAxis{10001, -5000.5, 5000.5, "PDG code"};

    AxisSpec dEdxAxis{400, 0., 200., "TPC dE/dx (arb. units)"};
    AxisSpec betaAxis{160, 0., 1.6, "#beta_{TOF}"};

    // Invariant-mass spectra
    histos.add("hInvMassPKUnlike",
               "pK invariant mass (unlike-sign);M_{pK} (GeV/c^{2});Counts",
               HistType::kTH1F, {massAxis});
    histos.add("hInvMassPKLike",
               "pK invariant mass (like-sign);M_{pK} (GeV/c^{2});Counts",
               HistType::kTH1F, {massAxis});

    // Invariant-mass vs pair pT (use p_{T} of pK system)
    histos.add("hInvMassPKUnlikeVsPt",
               "pK invariant mass vs p_{T} (unlike-sign);M_{pK} (GeV/c^{2});p_{T}^{pK} (GeV/c);Counts",
               HistType::kTH2F, {massAxis, ptAxis});
    histos.add("hInvMassPKLikeVsPt",
               "pK invariant mass vs p_{T} (like-sign);M_{pK} (GeV/c^{2});p_{T}^{pK} (GeV/c);Counts",
               HistType::kTH2F, {massAxis, ptAxis});

    // THnSparse for invariant-mass analysis (mass, pT, y, phi)
    if (lstarEnableSparse.value != 0) {
      histos.add("hLambdaStarPKUnlikeSparse",
                 "#Lambda^{*}(1520) pK unlike-sign candidates;M_{pK} (GeV/c^{2});p_{T}^{pK} (GeV/c);y_{pK};#varphi_{pK}",
                 HistType::kTHnSparseF,
                 {AxisSpec{400, 1.3, 1.9, "M_{pK} (GeV/c^{2})"},
                  AxisSpec{100, 0., 10., "p_{T}^{pK} (GeV/c)"},
                  AxisSpec{60, -1.5, 1.5, "y_{pK}"},
                  AxisSpec{64, -3.2, 3.2, "#varphi_{pK}"}, centAxis});

      histos.add("hLambdaStarPKLikeSparse",
                 "#Lambda^{*}(1520) pK like-sign candidates;M_{pK} (GeV/c^{2});p_{T}^{pK} (GeV/c);y_{pK};#varphi_{pK}",
                 HistType::kTHnSparseF,
                 {AxisSpec{400, 1.3, 1.9, "M_{pK} (GeV/c^{2})"},
                  AxisSpec{100, 0., 10., "p_{T}^{pK} (GeV/c)"},
                  AxisSpec{60, -1.5, 1.5, "y_{pK}"},
                  AxisSpec{64, -3.2, 3.2, "#varphi_{pK}"}, centAxis});

      histos.add("hLambdaStarPKMixedSparse",
                 "#Lambda^{*}(1520) pK mixed-event candidates;M_{pK} (GeV/c^{2});p_{T}^{pK} (GeV/c);y_{pK};#varphi_{pK}",
                 HistType::kTHnSparseF,
                 {AxisSpec{400, 1.3, 1.9, "M_{pK} (GeV/c^{2})"},
                  AxisSpec{100, 0., 10., "p_{T}^{pK} (GeV/c)"},
                  AxisSpec{60, -1.5, 1.5, "y_{pK}"},
                  AxisSpec{64, -3.2, 3.2, "#varphi_{pK}"}, centAxis});

      // THnSparse for deuteron-proxy invariant-mass analysis (mass, pT, y, phi)
      histos.add("hLambdaStarProxySparse",
                 "#Lambda^{*}(1520) deuteron-proxy candidates;M_{p_{proxy}K} (GeV/c^{2});p_{T}^{p_{proxy}K} (GeV/c);y_{p_{proxy}K};#varphi_{p_{proxy}K}",
                 HistType::kTHnSparseF,
                 {AxisSpec{400, 1.3, 1.9, "M_{p_{proxy}K} (GeV/c^{2})"},
                  AxisSpec{100, 0., 10., "p_{T}^{p_{proxy}K} (GeV/c)"},
                  AxisSpec{60, -1.5, 1.5, "y_{p_{proxy}K}"},
                  AxisSpec{64, -3.2, 3.2, "#varphi_{p_{proxy}K}"}, centAxis});

      histos.add("hLambdaStarProxyMixedSparse",
                 "#Lambda^{*}(1520) deuteron-proxy mixed-event candidates;M_{p_{proxy}K} (GeV/c^{2});p_{T}^{p_{proxy}K} (GeV/c);y_{p_{proxy}K};#varphi_{p_{proxy}K}",
                 HistType::kTHnSparseF,
                 {AxisSpec{400, 1.3, 1.9, "M_{p_{proxy}K} (GeV/c^{2})"},
                  AxisSpec{100, 0., 10., "p_{T}^{p_{proxy}K} (GeV/c)"},
                  AxisSpec{60, -1.5, 1.5, "y_{p_{proxy}K}"},
                  AxisSpec{64, -3.2, 3.2, "#varphi_{p_{proxy}K}"}, centAxis});
    }

    // Deuteron-proxy invariant mass (p_{proxy} from d/2 combined with K)
    histos.add("hDeuteronProxyMass",
               "#Lambda^{*} proxy invariant mass from (d/2 + K);M_{pK} (GeV/c^{2});Counts",
               HistType::kTH1F, {massAxis});

    // TPC dE/dx vs total momentum
    histos.add("hTPCdEdxVsP",
               "TPC dE/dx vs p;p (GeV/c);dE/dx (arb. units);Counts",
               HistType::kTH2F, {pAxis, dEdxAxis});

    // TOF #beta vs total momentum
    histos.add("hTOFBetaVsP",
               "TOF #beta vs p;p (GeV/c);#beta_{TOF};Counts",
               HistType::kTH2F, {pAxis, betaAxis});

    // --- Per-species PID QA (tagged) ---
    histos.add("hTPCdEdxVsP_Pi",
               "TPC dE/dx vs p (tagged #pi);p (GeV/c);dE/dx (arb. units);Counts",
               HistType::kTH2F, {pAxis, dEdxAxis});
    histos.add("hTPCdEdxVsP_K",
               "TPC dE/dx vs p (tagged K);p (GeV/c);dE/dx (arb. units);Counts",
               HistType::kTH2F, {pAxis, dEdxAxis});
    histos.add("hTPCdEdxVsP_P",
               "TPC dE/dx vs p (tagged p);p (GeV/c);dE/dx (arb. units);Counts",
               HistType::kTH2F, {pAxis, dEdxAxis});
    histos.add("hTPCdEdxVsP_D",
               "TPC dE/dx vs p (tagged d);p (GeV/c);dE/dx (arb. units);Counts",
               HistType::kTH2F, {pAxis, dEdxAxis});

    histos.add("hTOFBetaVsP_Pi",
               "TOF #beta vs p (tagged #pi);p (GeV/c);#beta_{TOF};Counts",
               HistType::kTH2F, {pAxis, betaAxis});
    histos.add("hTOFBetaVsP_K",
               "TOF #beta vs p (tagged K);p (GeV/c);#beta_{TOF};Counts",
               HistType::kTH2F, {pAxis, betaAxis});
    histos.add("hTOFBetaVsP_P",
               "TOF #beta vs p (tagged p);p (GeV/c);#beta_{TOF};Counts",
               HistType::kTH2F, {pAxis, betaAxis});
    histos.add("hTOFBetaVsP_D",
               "TOF #beta vs p (tagged d);p (GeV/c);#beta_{TOF};Counts",
               HistType::kTH2F, {pAxis, betaAxis});

    // --- MC QA (AO2D-MC mode) ---
    histos.add("hMcPrimariesPtEta",
               "MC charged physical primaries; p_{T} (GeV/c); #eta; Counts",
               HistType::kTH2F, {ptAxis, etaAxis});

    histos.add("hRecoMatchedPdg",
               "Reco tracks matched to MC (PDG); PDG code; Counts",
               HistType::kTH1F, {pdgAxis});

    histos.add("hRecoFakePtEta",
               "Reco tracks without MC label (fakes / unmatched); p_{T} (GeV/c); #eta; Counts",
               HistType::kTH2F, {ptAxis, etaAxis});

    // Deuteron-proxy kinematics and PID QA
    histos.add("hDeuteronProxyPt",
               "Deuteron proxy p_{T};p_{T} (GeV/c);Counts",
               HistType::kTH1F, {ptAxis});
    histos.add("hDeuteronProxyEta",
               "Deuteron proxy #eta;#eta;Counts",
               HistType::kTH1F, {etaAxis});
    histos.add("hDeuteronProxyPhi",
               "Deuteron proxy #varphi;#varphi;Counts",
               HistType::kTH1F, {phiAxis});

    histos.add("hNsigmaTPCDeuteron",
               "TPC n#sigma_{d};n#sigma^{TPC}_{d};Counts",
               HistType::kTH1F, {nsAxis});
    histos.add("hNsigmaTOFDeuteron",
               "TOF n#sigma_{d};n#sigma^{TOF}_{d};Counts",
               HistType::kTH1F, {nsAxis});

    histos.add("hNsigmaTPCDeuteronVsP",
               "TPC n#sigma_{d} vs p; p (GeV/c); n#sigma^{TPC}_{d};Counts",
               HistType::kTH2F, {pAxis, nsAxis});
    histos.add("hNsigmaTOFDeuteronVsP",
               "TOF n#sigma_{d} vs p; p (GeV/c); n#sigma^{TOF}_{d};Counts",
               HistType::kTH2F, {pAxis, nsAxis});

    // Kaon kinematics and PID QA
    histos.add("hKaonPt",
               "Kaon p_{T};p_{T} (GeV/c);Counts",
               HistType::kTH1F, {ptAxis});
    histos.add("hKaonEta",
               "Kaon #eta;#eta;Counts",
               HistType::kTH1F, {etaAxis});
    histos.add("hKaonPhi",
               "Kaon #varphi;#varphi;Counts",
               HistType::kTH1F, {phiAxis});

    histos.add("hNsigmaTPCKaon",
               "TPC n#sigma_{K};n#sigma^{TPC}_{K};Counts",
               HistType::kTH1F, {nsAxis});
    histos.add("hNsigmaTOFKaon",
               "TOF n#sigma_{K};n#sigma^{TOF}_{K};Counts",
               HistType::kTH1F, {nsAxis});

    histos.add("hNsigmaTPCKaonVsP",
               "TPC n#sigma_{K} vs p; p (GeV/c); n#sigma^{TPC}_{K};Counts",
               HistType::kTH2F, {pAxis, nsAxis});
    histos.add("hNsigmaTOFKaonVsP",
               "TOF n#sigma_{K} vs p; p (GeV/c); n#sigma^{TOF}_{K};Counts",
               HistType::kTH2F, {pAxis, nsAxis});
  }

  // AO2D-MC QA: truth primaries + reco-to-MC matching sanity plots
  void processMCQA(FilteredCollisionsMC::iterator const& collision,
                   FilteredTracksMC const& tracks,
                   aod::McParticles const& mcParticles)
  {
    if (lstarFillMCTruth.value == 0) {
      return;
    }

    // Same basic event selection as AO2D reco
    if (!keepCollisionAO2D(collision)) {
      return;
    }
    if (std::abs(collision.posZ()) > lstarCutVertex.value) {
      return;
    }

    // Require a linked MC collision
    if constexpr (requires { collision.has_mcCollision(); }) {
      if (!collision.has_mcCollision()) {
        return;
      }
    }

    int mcCollIdx = -1;
    if constexpr (requires { collision.mcCollisionId(); }) {
      mcCollIdx = collision.mcCollisionId();
    } else if constexpr (requires { collision.mcCollision().globalIndex(); }) {
      mcCollIdx = collision.mcCollision().globalIndex();
    }
    if (mcCollIdx < 0) {
      return;
    }

    // --- Truth QA: charged physical primaries ---
    auto truth = mcParticles.sliceBy(mcParticlesPerCollision, mcCollIdx);
    for (auto const& mcPart : truth) {
      if (!mcPart.isPhysicalPrimary()) {
        continue;
      }
      if (mcPart.pt() < lstarCutMCPtMin.value || std::abs(mcPart.eta()) > lstarCutMCEtaMax.value) {
        continue;
      }
      if (chargeFromPdg(mcPart.pdgCode()) == 0) {
        continue;
      }
      histos.fill(HIST("hMcPrimariesPtEta"), mcPart.pt(), mcPart.eta());
    }

    // --- Reco->MC matching QA ---
    const auto collIdx = collision.globalIndex();
    for (auto const& trk : tracks) {
      if (trk.collisionId() != collIdx) {
        continue;
      }
      if (!passTrackQuality(trk)) {
        continue;
      }

      if constexpr (requires { trk.has_mcParticle(); trk.mcParticle(); }) {
        if (trk.has_mcParticle()) {
          const auto mcPart = trk.mcParticle();
          histos.fill(HIST("hRecoMatchedPdg"), mcPart.pdgCode());
        } else {
          histos.fill(HIST("hRecoFakePtEta"), trk.pt(), trk.eta());
        }
      }
    }
  }

  PROCESS_SWITCH(Lambdastarproxy, processMCQA,
                 "AO2D-MC: fill truth and reco-matching QA histograms", false);

  // Helper: fill TPC dE/dx vs total momentum if TPC signal is available
  template <typename TTrack>
  void fillTPCdEdxVsPIfAvailable(const TTrack& trk)
  {
    if (lstarEnablePidQA.value == 0) {
      return;
    }
    // aod::TracksExtra provides tpcSignal(); keep the constexpr-guard for robustness
    if constexpr (requires { trk.tpcSignal(); }) {
      const float p = trk.p();
      histos.fill(HIST("hTPCdEdxVsP"), p, trk.tpcSignal());
    }
  }

  // Helper: fill TOF beta vs total momentum if beta is available (and looks valid)
  template <typename TTrack>
  void fillTOFBetaVsPIfAvailable(const TTrack& trk)
  {
    if (lstarEnablePidQA.value == 0) {
      return;
    }
    if constexpr (requires { trk.beta(); }) {
      bool hasTof = true;
      if constexpr (requires { trk.hasTOF(); }) {
        hasTof = trk.hasTOF();
      }
      const float beta = trk.beta();
      // Guard against default/invalid values for tracks without TOF match
      if (hasTof && beta > TofBetaMin && beta < TofBetaMax) {
        histos.fill(HIST("hTOFBetaVsP"), trk.p(), beta);
      }
    }
  }

  // --- Per-species PID QA helpers ---
  template <typename TTrack>
  static bool hasTOFMatch(const TTrack& trk)
  {
    if constexpr (requires { trk.hasTOF(); }) {
      return trk.hasTOF();
    }
    return true; // fallback: if column not present, assume available
  }

  // Return: 0=#pi, 1=K, 2=p, 3=d, -1=unclassified
  template <typename TTrack>
  int classifyPidSpecies(const TTrack& trk)
  {
    const bool hasTof = hasTOFMatch(trk);

    auto score = [hasTof](float nsTPC, float nsTOF) {
      return std::abs(nsTPC) + (hasTof ? std::abs(nsTOF) : 0.f);
    };

    float bestScore = 1e9f;
    int best = -1;

    // pion
    {
      const float nsTPC = trk.tpcNSigmaPi();
      const float nsTOF = trk.tofNSigmaPi();
      const bool pass = (std::abs(nsTPC) < lstarCutNsigmaTPCPi.value) && (!hasTof || (std::abs(nsTOF) < lstarCutNsigmaTOFPi.value));
      if (pass) {
        const float sc = score(nsTPC, nsTOF);
        if (sc < bestScore) {
          bestScore = sc;
          best = 0;
        }
      }
    }

    // kaon
    {
      const float nsTPC = trk.tpcNSigmaKa();
      const float nsTOF = trk.tofNSigmaKa();
      const bool pass = (std::abs(nsTPC) < lstarCutNsigmaTPCKaon.value) && (!hasTof || (std::abs(nsTOF) < lstarCutNsigmaTOFKaon.value));
      if (pass) {
        const float sc = score(nsTPC, nsTOF);
        if (sc < bestScore) {
          bestScore = sc;
          best = 1;
        }
      }
    }

    // proton
    {
      const float nsTPC = trk.tpcNSigmaPr();
      const float nsTOF = trk.tofNSigmaPr();
      const bool pass = (std::abs(nsTPC) < lstarCutNsigmaTPCPr.value) && (!hasTof || (std::abs(nsTOF) < lstarCutNsigmaTOFPr.value));
      if (pass) {
        const float sc = score(nsTPC, nsTOF);
        if (sc < bestScore) {
          bestScore = sc;
          best = 2;
        }
      }
    }

    // deuteron
    {
      const float nsTPC = trk.tpcNSigmaDe();
      const float nsTOF = trk.tofNSigmaDe();
      const bool pass = (std::abs(nsTPC) < lstarCutNsigmaTPCDe.value) && (!hasTof || (std::abs(nsTOF) < lstarCutNsigmaTOFDe.value));
      if (pass) {
        const float sc = score(nsTPC, nsTOF);
        if (sc < bestScore) {
          bestScore = sc;
          best = 3;
        }
      }
    }

    return best;
  }

  // Helper to compute invariant mass from two 3-momenta and masses
  static double invariantMass(float px1, float py1, float pz1, double m1,
                              float px2, float py2, float pz2, double m2)
  {
    const double e1 = std::sqrt(m1 * m1 + px1 * px1 + py1 * py1 + pz1 * pz1);
    const double e2 = std::sqrt(m2 * m2 + px2 * px2 + py2 * py2 + pz2 * pz2);
    const double ex = px1 + px2;
    const double ey = py1 + py2;
    const double ez = pz1 + pz2;
    const double eTot = e1 + e2;
    const double p2Tot = ex * ex + ey * ey + ez * ez;
    const double m2Tot = eTot * eTot - p2Tot;
    return m2Tot > 0. ? std::sqrt(m2Tot) : 0.;
  }

  void process(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks)
  {
    // Event selection (cfgTrigger) -- AO2D only
    if (!keepCollisionAO2D(collision)) {
      return;
    }
    // physics masses (GeV/c^2)
    constexpr double MassProton = o2::constants::physics::MassProton;
    constexpr double MassKaonCharged = o2::constants::physics::MassKaonCharged;

    // --- Inclusive PID QA: keep #pi/K/p/d bands in TPC dE/dx and TOF beta plots ---
    if (lstarEnablePidQA.value != 0) {
      for (auto const& trk : tracks) {
        if (trk.pt() < lstarCutPtMin.value || std::abs(trk.eta()) > lstarCutEtaMax.value) {
          continue;
        }
        if (!passTrackQuality(trk)) {
          continue;
        }
        if (trk.sign() == 0) {
          continue;
        }
        // Inclusive PID QA
        fillTPCdEdxVsPIfAvailable(trk);
        fillTOFBetaVsPIfAvailable(trk);

        // Per-species PID-QA (tagged) histograms
        const int sp = classifyPidSpecies(trk);
        switch (sp) {
          case 0: { // pion
            if constexpr (requires { trk.tpcSignal(); }) {
              histos.fill(HIST("hTPCdEdxVsP_Pi"), trk.p(), trk.tpcSignal());
            }
            if constexpr (requires { trk.beta(); }) {
              const bool hasTof = hasTOFMatch(trk);
              const float beta = trk.beta();
              if (hasTof && beta > TofBetaMin && beta < TofBetaMax) {
                histos.fill(HIST("hTOFBetaVsP_Pi"), trk.p(), beta);
              }
            }
            break;
          }

          case 1: { // kaon
            if constexpr (requires { trk.tpcSignal(); }) {
              histos.fill(HIST("hTPCdEdxVsP_K"), trk.p(), trk.tpcSignal());
            }
            if constexpr (requires { trk.beta(); }) {
              const bool hasTof = hasTOFMatch(trk);
              const float beta = trk.beta();
              if (hasTof && beta > TofBetaMin && beta < TofBetaMax) {
                histos.fill(HIST("hTOFBetaVsP_K"), trk.p(), beta);
              }
            }
            break;
          }

          case 2: { // proton
            if constexpr (requires { trk.tpcSignal(); }) {
              histos.fill(HIST("hTPCdEdxVsP_P"), trk.p(), trk.tpcSignal());
            }
            if constexpr (requires { trk.beta(); }) {
              const bool hasTof = hasTOFMatch(trk);
              const float beta = trk.beta();
              if (hasTof && beta > TofBetaMin && beta < TofBetaMax) {
                histos.fill(HIST("hTOFBetaVsP_P"), trk.p(), beta);
              }
            }
            break;
          }

          case 3: { // deuteron
            if constexpr (requires { trk.tpcSignal(); }) {
              histos.fill(HIST("hTPCdEdxVsP_D"), trk.p(), trk.tpcSignal());
            }
            if constexpr (requires { trk.beta(); }) {
              const bool hasTof = hasTOFMatch(trk);
              const float beta = trk.beta();
              if (hasTof && beta > TofBetaMin && beta < TofBetaMax) {
                histos.fill(HIST("hTOFBetaVsP_D"), trk.p(), beta);
              }
            }
            break;
          }

          default:
            break;
        }
      }
    }

    std::vector<KaonCand> kaonCands;
    std::vector<ProxyCand> proxyCands;
    kaonCands.reserve(128);
    proxyCands.reserve(32);

    float eventMultFallback = 0.f; // fallback mixing variable: number of selected charged tracks (after quality cuts)

    // Inclusive PID QA loop: count all selected charged tracks for fallback multiplicity
    for (auto const& trk : tracks) {
      if (trk.pt() < lstarCutPtMin.value || std::abs(trk.eta()) > lstarCutEtaMax.value) {
        continue;
      }
      if (!passTrackQuality(trk)) {
        continue;
      }
      if (trk.sign() == 0) {
        continue;
      }
      eventMultFallback += 1.f;
      if (lstarEnablePidQA.value == 0) {
        continue;
      }
      // Inclusive PID QA
      fillTPCdEdxVsPIfAvailable(trk);
      fillTOFBetaVsPIfAvailable(trk);

      // Per-species PID-QA (tagged) histograms
      const int sp = classifyPidSpecies(trk);
      switch (sp) {
        case 0: { // pion
          if constexpr (requires { trk.tpcSignal(); }) {
            histos.fill(HIST("hTPCdEdxVsP_Pi"), trk.p(), trk.tpcSignal());
          }
          if constexpr (requires { trk.beta(); }) {
            const bool hasTof = hasTOFMatch(trk);
            const float beta = trk.beta();
            if (hasTof && beta > TofBetaMin && beta < TofBetaMax) {
              histos.fill(HIST("hTOFBetaVsP_Pi"), trk.p(), beta);
            }
          }
          break;
        }
        case 1: { // kaon
          if constexpr (requires { trk.tpcSignal(); }) {
            histos.fill(HIST("hTPCdEdxVsP_K"), trk.p(), trk.tpcSignal());
          }
          if constexpr (requires { trk.beta(); }) {
            const bool hasTof = hasTOFMatch(trk);
            const float beta = trk.beta();
            if (hasTof && beta > TofBetaMin && beta < TofBetaMax) {
              histos.fill(HIST("hTOFBetaVsP_K"), trk.p(), beta);
            }
          }
          break;
        }
        case 2: { // proton
          if constexpr (requires { trk.tpcSignal(); }) {
            histos.fill(HIST("hTPCdEdxVsP_P"), trk.p(), trk.tpcSignal());
          }
          if constexpr (requires { trk.beta(); }) {
            const bool hasTof = hasTOFMatch(trk);
            const float beta = trk.beta();
            if (hasTof && beta > TofBetaMin && beta < TofBetaMax) {
              histos.fill(HIST("hTOFBetaVsP_P"), trk.p(), beta);
            }
          }
          break;
        }
        case 3: { // deuteron
          if constexpr (requires { trk.tpcSignal(); }) {
            histos.fill(HIST("hTPCdEdxVsP_D"), trk.p(), trk.tpcSignal());
          }
          if constexpr (requires { trk.beta(); }) {
            const bool hasTof = hasTOFMatch(trk);
            const float beta = trk.beta();
            if (hasTof && beta > TofBetaMin && beta < TofBetaMax) {
              histos.fill(HIST("hTOFBetaVsP_D"), trk.p(), beta);
            }
          }
          break;
        }
        default:
          break;
      }
    }

    // Compute event multiplicity (FT0M or fallback)
    const float eventMult = eventMultiplicityFT0MOrFallback(collision, eventMultFallback);

    // Deuteron candidates -> proton-proxy candidates
    for (auto const& trkD : tracks) {
      if (trkD.pt() < lstarCutPtMin.value || std::abs(trkD.eta()) > lstarCutEtaMax.value) {
        continue;
      }
      if (!passTrackQuality(trkD)) {
        continue;
      }
      if (trkD.sign() == 0) {
        continue;
      }

      // PID for deuteron candidates
      const float nsTPCDe = trkD.tpcNSigmaDe();
      const float nsTOFDe = trkD.tofNSigmaDe();
      const bool isDeuteron = (std::abs(nsTPCDe) < lstarCutNsigmaTPCDe.value) && (std::abs(nsTOFDe) < lstarCutNsigmaTOFDe.value);
      if (!isDeuteron) {
        continue;
      }

      // Deuteron kinematics
      const float ptD = trkD.pt();
      const float etaD = trkD.eta();
      const float phiD = trkD.phi();
      const double pD = static_cast<double>(ptD) * std::cosh(static_cast<double>(etaD));

      // QA histos for deuteron PID and kinematics
      if (lstarEnablePidQA.value != 0) {
        histos.fill(HIST("hDeuteronProxyPt"), ptD);
        histos.fill(HIST("hDeuteronProxyEta"), etaD);
        histos.fill(HIST("hDeuteronProxyPhi"), phiD);
        histos.fill(HIST("hNsigmaTPCDeuteron"), nsTPCDe);
        histos.fill(HIST("hNsigmaTOFDeuteron"), nsTOFDe);
        histos.fill(HIST("hNsigmaTPCDeuteronVsP"), pD, nsTPCDe);
        histos.fill(HIST("hNsigmaTOFDeuteronVsP"), pD, nsTOFDe);
      }

      // build proton-proxy momentum from deuteron: p_p ≈ p_d / 2
      const float pxProxy = ProxyMomentumScale * ptD * std::cos(phiD);
      const float pyProxy = ProxyMomentumScale * ptD * std::sin(phiD);
      const float pzProxy = ProxyMomentumScale * ptD * std::sinh(etaD);

      proxyCands.push_back(ProxyCand{pxProxy, pyProxy, pzProxy, static_cast<int>(trkD.sign()), static_cast<int>(trkD.globalIndex())});
    }

    // Kaon candidates
    for (auto const& trkK : tracks) {
      if (trkK.pt() < lstarCutPtMin.value || std::abs(trkK.eta()) > lstarCutEtaMax.value) {
        continue;
      }
      if (!passTrackQuality(trkK)) {
        continue;
      }
      if (trkK.sign() == 0) {
        continue;
      }

      // PID for kaon candidates
      const float nsTPCK = trkK.tpcNSigmaKa();
      const float nsTOFK = trkK.tofNSigmaKa();
      const bool isKaon = (std::abs(nsTPCK) < lstarCutNsigmaTPCKaon.value) && (std::abs(nsTOFK) < lstarCutNsigmaTOFKaon.value);
      if (!isKaon) {
        continue;
      }

      // Kaon kinematics
      const float ptK = trkK.pt();
      const float etaK = trkK.eta();
      const float phiK = trkK.phi();
      const double pK = static_cast<double>(ptK) * std::cosh(static_cast<double>(etaK));

      // Kaon QA
      if (lstarEnablePidQA.value != 0) {
        histos.fill(HIST("hKaonPt"), ptK);
        histos.fill(HIST("hKaonEta"), etaK);
        histos.fill(HIST("hKaonPhi"), phiK);
        histos.fill(HIST("hNsigmaTPCKaon"), nsTPCK);
        histos.fill(HIST("hNsigmaTOFKaon"), nsTOFK);
        histos.fill(HIST("hNsigmaTPCKaonVsP"), pK, nsTPCK);
        histos.fill(HIST("hNsigmaTOFKaonVsP"), pK, nsTOFK);
      }

      const float pxK = ptK * std::cos(phiK);
      const float pyK = ptK * std::sin(phiK);
      const float pzK = ptK * std::sinh(etaK);

      kaonCands.push_back(KaonCand{pxK, pyK, pzK, static_cast<int>(trkK.sign()), static_cast<int>(trkK.globalIndex())});
    }

    if (proxyCands.empty() || kaonCands.empty()) {
      // still update mixing buffer so that later events can mix with this one
      LStarMixEventEntry entry;
      entry.mult = eventMult;
      entry.zvtx = collision.posZ();
      entry.kaons = std::move(kaonCands);
      entry.proxies = std::move(proxyCands);
      mLStarMixEvents.push_front(std::move(entry));
      if (mLStarMixEvents.size() > static_cast<size_t>(lstarNoMixedEvents.value)) {
        mLStarMixEvents.pop_back();
      }
      return;
    }

    // --- SAME-EVENT: proxy (d/2) + K ---
    for (auto const& pr : proxyCands) {
      for (auto const& k : kaonCands) {
        if (pr.tid == k.tid)
          continue; // sanity check: should never match, but just in case of bug in candidate-building logic
        const double mass = invariantMass(pr.px, pr.py, pr.pz, MassProton, k.px, k.py, k.pz, MassKaonCharged);

        const float pxTot = pr.px + k.px;
        const float pyTot = pr.py + k.py;
        const float pzTot = pr.pz + k.pz;
        const float ptPair = std::sqrt(pxTot * pxTot + pyTot * pyTot);
        const float phiPair = phiFromPxPy(pxTot, pyTot);

        const double eTot = std::sqrt(mass * mass + static_cast<double>(pxTot) * pxTot + static_cast<double>(pyTot) * pyTot + static_cast<double>(pzTot) * pzTot);
        const float yPair = rapidityFromEPz(eTot, pzTot);

        // Inclusive invariant-mass spectrum for the #Lambda^{*} proxy
        histos.fill(HIST("hDeuteronProxyMass"), mass);
        if (lstarEnableSparse.value != 0) {
          histos.fill(HIST("hLambdaStarProxySparse"), mass, ptPair, yPair, phiPair, eventMult);
        }

        const bool unlikeSign = (pr.charge * k.charge) < 0;
        if (unlikeSign) {
          histos.fill(HIST("hInvMassPKUnlike"), mass);
          histos.fill(HIST("hInvMassPKUnlikeVsPt"), mass, ptPair);
          if (lstarEnableSparse.value != 0) {
            histos.fill(HIST("hLambdaStarPKUnlikeSparse"), mass, ptPair, yPair, phiPair, eventMult);
          }
        } else {
          histos.fill(HIST("hInvMassPKLike"), mass);
          histos.fill(HIST("hInvMassPKLikeVsPt"), mass, ptPair);
          if (lstarEnableSparse.value != 0) {
            histos.fill(HIST("hLambdaStarPKLikeSparse"), mass, ptPair, yPair, phiPair, eventMult);
          }
        }
      }
    }

    // --- MIXED-EVENT: current proxies + previous-event kaons ---
    for (auto const& prev : mLStarMixEvents) {
      if (std::abs(prev.zvtx - collision.posZ()) > lstarMixZvtxMax.value)
        continue;
      if (std::abs(prev.mult - eventMult) > lstarMixMultMax.value)
        continue;
      if (prev.kaons.empty()) {
        continue;
      }

      for (auto const& pr : proxyCands) {
        for (auto const& k : prev.kaons) {
          // convention: mix for unlike-sign only (resonance background)
          if ((pr.charge * k.charge) >= 0) {
            continue;
          }
          if (pr.tid == k.tid)
            continue; // sanity check: should never match, but just in case of bug in candidate-building logic

          const double mass = invariantMass(pr.px, pr.py, pr.pz, MassProton, k.px, k.py, k.pz, MassKaonCharged);

          const float pxTot = pr.px + k.px;
          const float pyTot = pr.py + k.py;
          const float pzTot = pr.pz + k.pz;
          const float ptPair = std::sqrt(pxTot * pxTot + pyTot * pyTot);
          const float phiPair = phiFromPxPy(pxTot, pyTot);

          const double eTot = std::sqrt(mass * mass + static_cast<double>(pxTot) * pxTot + static_cast<double>(pyTot) * pyTot + static_cast<double>(pzTot) * pzTot);
          const float yPair = rapidityFromEPz(eTot, pzTot);

          // Fill mixed-event THnSparse
          if (lstarEnableSparse.value != 0) {
            histos.fill(HIST("hLambdaStarPKMixedSparse"), mass, ptPair, yPair, phiPair, eventMult);
            histos.fill(HIST("hLambdaStarProxyMixedSparse"), mass, ptPair, yPair, phiPair, eventMult);
          }
        }
      }
    }

    // --- Update mixing buffer with current event ---
    LStarMixEventEntry entry;
    entry.mult = eventMult;
    entry.zvtx = collision.posZ();
    entry.kaons = std::move(kaonCands);
    entry.proxies = std::move(proxyCands);

    mLStarMixEvents.push_front(std::move(entry));
    if (mLStarMixEvents.size() > static_cast<size_t>(lstarNoMixedEvents.value)) {
      mLStarMixEvents.pop_back();
    }
  }

  PROCESS_SWITCH(Lambdastarproxy, process, "Lambda* proxy via (d/2)+K", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Nucleibalance>(cfgc),
    adaptAnalysisTask<Lambdastarproxy>(cfgc)};
}
