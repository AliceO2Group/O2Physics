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

/// \file lambdastarproxy.cxx
/// \brief Standalone Lambda* deuteron-proxy task extracted from nucleibalance.cxx.
/// \author Sushanta Tripathy <sushanta.tripathy@cern.ch>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <deque>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

struct Lambdastarproxy {
  Service<o2::framework::O2DatabasePDG> pdgDB{};

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

  static constexpr int FT0MPercentileBins = 100;
  static constexpr float FT0MPercentileMax = 100.f;
  enum class FT0MEventStatus : int {
    Valid = 1,
    Negative,
    OutOfRange,
    NonFinite,
    Unavailable
  };

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
  // PID strategy values
  static constexpr int PidStrategyRectangular = 0;
  static constexpr int PidStrategyCircularTPCAndTOF = 1;
  static constexpr int PidStrategyNucleiDeuteronTPC = 2;
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
  Configurable<int> lstarEnableTOFNsigmaCutDe{"lstarEnableTOFNsigmaCutDe", 0, "Enable deuteron-only TOF nSigma cut in PID strategy 2"};
  Configurable<int> lstarProxyUseProtonPIDAsDeuteron{"lstarProxyUseProtonPIDAsDeuteron", 0, "Closure/control test: keep proxy candidates that pass the normal deuteron selection and are also proton-like, then build p_proxy = p_track/2"};
  // Optional deuteron-only TOF auxiliary selections.
  // Defaults are OFF, so strategy 2 remains TPC nSigma only for deuterons.
  Configurable<int> lstarEnableBetaCutDe{"lstarEnableBetaCutDe", 0, "Enable deuteron-only TOF beta cut using beta() > lstarBetaCutDe"};
  Configurable<float> lstarBetaCutDe{"lstarBetaCutDe", 0.4f, "Minimum TOF beta for deuteron-only beta cut"};
  Configurable<int> lstarEnableExpSignalTOFDe{"lstarEnableExpSignalTOFDe", 0, "Enable deuteron-only TOF expected-signal-difference cut when lstarTOFExpSignalDiffDeMax > 0"};
  Configurable<float> lstarTOFExpSignalDiffDeMax{"lstarTOFExpSignalDiffDeMax", -1.f, "Maximum |tofExpSignalDiffDe| for deuterons; <=0 disables the numeric cut"};

  // PID strategy for final K/p/d candidate selection.
  // 0 = pT-ref dependent rectangular cuts:
  //     pT < pTref: require |TPC| < TPC cut only
  //     pT >= pTref and TOF exists: require |TPC| < TPC cut and |TOF| < TOF cut
  //     pT >= pTref and TOF missing: reject
  // 1 = pT-ref dependent circular cut:
  //     pT < pTref: require |TPC| < TPC cut only
  //     pT >= pTref and TOF exists: require sqrt(TPC^2 + TOF^2) < circular cut
  //     pT >= pTref and TOF missing: reject
  // 2 = hybrid official-like PID:
  //     deuterons: require |TPC_d| < TPC cut only
  //     kaons/protons: use the Lambda(1520)-like pT-ref circular TPC+TOF logic
  Configurable<int> lstarPidStrategy{"lstarPidStrategy", int{PidStrategyRectangular}, "PID strategy: 0=pTref rectangular TPC/TOF, 1=pTref circular TPC+TOF, 2=hybrid: K/p circular, d TPC-only"};

  Configurable<float> lstarPidCircularCutKaon{"lstarPidCircularCutKaon", 2.0f, "Circular PID cut sqrt(nSigmaTPC_K^2+nSigmaTOF_K^2) for kaons"};
  Configurable<float> lstarPidCircularCutPr{"lstarPidCircularCutPr", 2.0f, "Circular PID cut sqrt(nSigmaTPC_p^2+nSigmaTOF_p^2) for protons"};
  Configurable<float> lstarPidCircularCutDe{"lstarPidCircularCutDe", 2.0f, "Circular PID cut sqrt(nSigmaTPC_d^2+nSigmaTOF_d^2) for deuterons"};

  Configurable<float> lstarPidPtRefKaon{"lstarPidPtRefKaon", 0.5f, "pT reference for kaon PID strategy"};
  Configurable<float> lstarPidPtRefPr{"lstarPidPtRefPr", 0.8f, "pT reference for proton PID strategy"};
  Configurable<float> lstarPidPtRefDe{"lstarPidPtRefDe", 0.8f, "pT reference for deuteron PID strategy"};

  // Track quality
  Configurable<int> lstarRequireINELgt0{"lstarRequireINELgt0", 1, "Require INEL>0 event selection using isInelGt0() when available; fallback to multNTracksPVeta1()/numContrib()"};
  Configurable<bool> lstarRequireGlobalTrack{"lstarRequireGlobalTrack", bool{RequireGlobalTrackDefault}, "Require global tracks (default)"};
  Configurable<int> lstarRequirePrimaryTrack{"lstarRequirePrimaryTrack", 1, "Require isPrimaryTrack() for Lambda* proxy candidates when the column is available"};
  Configurable<int> lstarOnlyGlobalTrackCuts{"lstarOnlyGlobalTrackCuts", 0, "If enabled, apply only the global-track plus primary-track requirements and skip additional ITS/TPC/DCA/chi2 track-quality cuts"};
  Configurable<int> lstarTPCNClsMin{"lstarTPCNClsMin", int{TPCNClsMinDefault}, "Minimum number of TPC crossed rows (tpcNClsCrossedRows)"};
  Configurable<float> lstarTPCCrossedRowsOverFindableMin{"lstarTPCCrossedRowsOverFindableMin", 0.8f, "Minimum TPC crossed rows over findable clusters"};
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
  // Master switch for PID QA histogram filling.
  // Inclusive PID-QA histograms are filled before final PID cuts, after event/track-quality cuts.
  // Candidate-level nSigma/TOF/DCA histograms are filled after final PID cuts.
  Configurable<int> lstarEnablePidQA{"lstarEnablePidQA", 0, "Enable PID QA histograms (dE/dx, TOF #beta, proxy invariant-mass QA, etc.): 1 = ON, 0 = OFF"};
  Configurable<int> lstarEnableProxyControls{"lstarEnableProxyControls", 1, "Enable TOF/PID-partitioned proxy spectra and inclusive full-pK maps (requires lstarEnableSparse)"};
  Configurable<int> lstarRequireTOFMatchDe{"lstarRequireTOFMatchDe", 0, "Require a TOF match for every selected deuteron; 0 preserves hybrid unmatched candidates"};
  Configurable<int> lstarEnableSparse{"lstarEnableSparse", 1, "Enable THnSparse invariant-mass histograms (#Lambda^{*} pK and proxy); 1 = ON, 0 = OFF"};
  Configurable<float> lstarLambdaAbsYMax{"lstarLambdaAbsYMax", 0.5f, "Max |y_{pK}| (or y_{proxy K}) for #Lambda^{*} candidates"};

  // Booking and filling use the same switches. Full mass maps grow with
  // occupied mass x mass x pT x multiplicity cells and are diagnostic only.
  Configurable<int> lstarEnableFullPKMaps{"lstarEnableFullPKMaps", 0, "Book/fill 4D and 8D proxy/full-pK maps; enable only for diagnostic productions"};
  Configurable<int> lstarEnableProtonLikeSparse{"lstarEnableProtonLikeSparse", 0, "Book/fill redundant proton-like proxy spectra (also obtainable from compact controls)"};
  Configurable<int> lstarSparseMultBins{"lstarSparseMultBins", 100, "Multiplicity-axis bins in all invariant-mass sparse histograms; keep identical across merge inputs"};
  Configurable<int> lstarMCPrimaryOnly{"lstarMCPrimaryOnly", 1, "Generated p/K/d fiducial spectra: require physical primaries (strong-decay daughters included)"};

  // Optional diagnostics leave the original selections and histograms intact.
  Configurable<int> lstarEnableMCProxySpecies{"lstarEnableMCProxySpecies", 0, "Book selected MC proxy species x TOF/proton/multiplicity states and retained-parent controls"};
  Configurable<int> lstarEnableProxyPIDProbes{"lstarEnableProxyPIDProbes", 0, "Book separate loose-TPC/TOF matched pair probes for PID-versus-proxy-mass fits (data or reconstructed MC)"};
  Configurable<float> lstarProxyProbeTPCDeMax{"lstarProxyProbeTPCDeMax", 6.f, "Loose d-slot TPC envelope for PID probes; 0 < cut <= 6"};
  Configurable<float> lstarProxyProbeTPCKaMax{"lstarProxyProbeTPCKaMax", 6.f, "Loose K-slot TPC envelope for PID probes; 0 < cut <= 6"};
  Configurable<float> lstarProxyProbeTOFDeMax{"lstarProxyProbeTOFDeMax", 8.f, "Loose d-slot TOF envelope for PID probes; 0 < cut <= 8"};
  Configurable<float> lstarProxyProbeTOFKaMax{"lstarProxyProbeTOFKaMax", 8.f, "Loose K-slot TOF envelope for PID probes; 0 < cut <= 8"};

  struct RetainedMother {
    int64_t id = -1;
    int pdg = 0;
  };

  struct TruthInfo {
    int64_t id = -1;
    int collisionId = -1;
    int pdg = 0;
    float px = 0.f, py = 0.f, pz = 0.f;
    bool primary = false;
    std::vector<int64_t> lambdaMothers;
    std::vector<RetainedMother> diagnosticMothers;
  };

  struct KaonCand {
    float px = 0.f, py = 0.f, pz = 0.f;
    bool hasTOF = false;
    int charge = 0;
    int tid = -1;
    TruthInfo truth{};
  };
  // PID flags describe detector compatibility, not truth-level purity.
  struct ProxyCand {
    float px = 0.f, py = 0.f, pz = 0.f;
    float pxFull = 0.f, pyFull = 0.f, pzFull = 0.f;
    bool protonLike = false;
    bool hasTOF = false;
    bool protonTPCCompatible = false;
    int charge = 0;
    int tid = -1;
    TruthInfo truth{};
  };
  struct ProtonCand {
    float px = 0.f, py = 0.f, pz = 0.f;
    int charge = 0;
    int tid = -1;
    TruthInfo truth{};
  };

  // Helpers for invariant-mass kinematics
  static float phiFromPxPy(float px, float py)
  {
    return RecoDecay::constrainAngle(std::atan2(py, px), -o2::constants::math::PI);
  }

  static float ptFromPxPy(float px, float py)
  {
    return RecoDecay::pt(std::array{px, py});
  }

  static float rapidityFromMomentumAndMass(float px, float py, float pz, double mass)

  {
    return RecoDecay::y(std::array{px, py, pz}, mass);
  }

  // Mixed-event pool entry for pK / proxy background (AO2D only)
  struct LStarMixEventEntry {
    float mult = 0.f;
    float zvtx = 0.f;
    std::vector<KaonCand> kaons;
    std::vector<ProtonCand> protons;
  };

  // Keep last N events for event-mixing
  std::deque<LStarMixEventEntry> mLStarMixEvents;

  template <typename TCollision>
  bool keepCollisionAO2D(TCollision const& collision) const
  {
    // Prefer the official event-selection flag when available.
    if (lstarRequireINELgt0.value != 0) {
      bool isINELgt0 = true;

      if constexpr (requires { collision.isInelGt0(); }) {
        isINELgt0 = collision.isInelGt0();
      } else if constexpr (requires { collision.multNTracksPVeta1(); }) {
        isINELgt0 = collision.multNTracksPVeta1() > 0;
      } else if constexpr (requires { collision.numContrib(); }) {
        isINELgt0 = collision.numContrib() > 0;
      }

      if (!isINELgt0) {
        return false;
      }
    }
    if (lstarCfgTrigger.value == TriggerNone) {
      return true;
    }
    if (lstarCfgTrigger.value == TriggerSel8) {
      return collision.sel8();
    }
    if (lstarCfgTrigger.value == TriggerSel8Quality) {
      return collision.sel8() &&
             collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
             collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
             collision.selection_bit(aod::evsel::kIsGoodITSLayersAll);
    }
    if (lstarCfgTrigger.value == TriggerSel8OccQuality) {
      const int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < lstarMinOcc.value || occupancy >= lstarMaxOcc.value) {
        return false;
      }
      return collision.sel8() &&
             collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
             collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
             collision.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) &&
             collision.selection_bit(aod::evsel::kIsGoodITSLayersAll);
    }
    if (lstarCfgTrigger.value == TriggerSel8NoSbpZvtx) {
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

  template <typename TCollision>
  void fillEventMultiplicityCounters(TCollision const& collision, float eventMult)
  {
    // Use the exact estimator of the pair maps, including the legacy track-count
    // fallback. Nonfinite values have no percentile and are counted in overflow.
    histos.fill(HIST("hEventMultiplicity"), std::isfinite(eventMult) ? eventMult : FT0MPercentileMax);

    // Keep the raw FT0M denominator separate: fallback track counts must never
    // acquire an FT0M percentile label. Missing/nonfinite FT0M goes to overflow;
    // finite negative and >=100 values retain their underflow/overflow meaning.
    float ft0m = FT0MPercentileMax;
    auto status = FT0MEventStatus::Unavailable;
    if constexpr (requires { collision.centFT0M(); }) {
      const float value = collision.centFT0M();
      if (!std::isfinite(value)) {
        status = FT0MEventStatus::NonFinite;
      } else {
        ft0m = value;
        if (value < 0.f) {
          status = FT0MEventStatus::Negative;
        } else if (value >= FT0MPercentileMax) {
          status = FT0MEventStatus::OutOfRange;
        } else {
          status = FT0MEventStatus::Valid;
        }
      }
    }
    histos.fill(HIST("hEventFT0M"), ft0m);
    histos.fill(HIST("hEventFT0MStatus"), static_cast<int>(status));
  }

  int chargeFromPdg(int pdg) const
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
        return (pdg >= 0) ? -1 : 1;
      case PdgPion:   // pi
      case PdgKaon:   // K
      case PdgProton: // p
        return (pdg >= 0) ? 1 : -1;
      default:
        if (auto* particle = pdgDB->GetParticle(pdg)) {
          return static_cast<int>(std::lround(particle->Charge() / 3.));
        }
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

    // Always require primary tracks
    if (lstarRequirePrimaryTrack.value != 0) {
      if constexpr (requires { trk.isPrimaryTrack(); }) {
        if (!trk.isPrimaryTrack()) {
          return false;
        }
      }
    }

    // Optional official-baseline mode: use only the O2 global-track and primary-track definitions.
    // This bypasses the extra explicit ITS/TPC/DCA/chi2 cuts below.
    if (lstarOnlyGlobalTrackCuts.value != 0) {
      return true;
    }

    if constexpr (requires { trk.itsNCls(); }) {
      if (lstarITSNClusters.value > 0 && trk.itsNCls() < lstarITSNClusters.value) {
        return false;
      }
    }

    if constexpr (requires { trk.tpcNClsCrossedRows(); }) {
      if (lstarTPCNClsMin.value > 0 && trk.tpcNClsCrossedRows() < lstarTPCNClsMin.value) {
        return false;
      }
    } else if constexpr (requires { trk.tpcNClsFindable(); trk.tpcCrossedRowsOverFindableCls(); }) {
      if (lstarTPCNClsMin.value > 0 && trk.tpcNClsFindable() * trk.tpcCrossedRowsOverFindableCls() < lstarTPCNClsMin.value) {
        return false;
      }
    }

    if constexpr (requires { trk.tpcCrossedRowsOverFindableCls(); }) {
      if (lstarTPCCrossedRowsOverFindableMin.value > 0.f &&
          trk.tpcCrossedRowsOverFindableCls() < lstarTPCCrossedRowsOverFindableMin.value) {
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
    if (doprocess && doprocessMCReco) {
      LOGF(fatal, "Enable either process (data) or processMCReco, never both");
    }
    if (lstarNoMixedEvents.value < 0 || lstarSparseMultBins.value <= 0) {
      LOGF(fatal, "lstarNoMixedEvents must be nonnegative; lstarSparseMultBins must be positive");
    }
    AxisSpec massAxis{200, 1.4, 1.9, "M_{pK} (GeV/c^{2})"};
    AxisSpec ptAxis{100, 0., 10., "p_{T} (GeV/c)"};
    AxisSpec nsAxis{500, -50., 50., "n#sigma"};
    AxisSpec tofMatchAxis{2, -0.5, 1.5, "has TOF match"};
    AxisSpec pAxis{100, 0., 10., "p (GeV/c)"};
    AxisSpec etaAxis{80, -2., 2., "#eta"};
    AxisSpec phiAxis{64, 0., o2::constants::math::TwoPI, "#varphi"};
    AxisSpec centAxis{lstarSparseMultBins.value, 0., 100., "multiplicity/centrality"};

    AxisSpec pdgAxis{10001, -5000.5, 5000.5, "PDG code"};

    AxisSpec dEdxAxis{400, 0., 200., "TPC dE/dx (arb. units)"};
    AxisSpec betaAxis{160, 0., 1.6, "#beta_{TOF}"};
    AxisSpec dcaXYAxis{200, -0.2, 0.2, "DCA_{xy} (cm)"};
    AxisSpec dcaZAxis{200, -0.2, 0.2, "DCA_{z} (cm)"};

    // Count reconstructed collisions received by the data or reco-MC process, before any
    // task event cuts and after all of them, independently of track candidates.
    histos.add("hEventSelection", "Proxy event selection;Selection stage;Events",
               HistType::kTH1D, {AxisSpec{2, 0.5, 2.5, "Selection stage"}});
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "Before event selection");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "After event selection");

    // Small event-level denominators, independent of candidate and QA switches.
    // Double storage and Sumw2 preserve large merged event counts and errors.
    histos.add("hEventFT0M", "Selected reconstructed events;Raw FT0M percentile;Events",
               HistType::kTH1D, {AxisSpec{FT0MPercentileBins, 0., FT0MPercentileMax, "Raw FT0M percentile"}}, true);
    histos.add("hEventMultiplicity", "Selected reconstructed events (same estimator as pair maps);Multiplicity/centrality;Events",
               HistType::kTH1D, {centAxis}, true);
    histos.add("hEventFT0MStatus", "Selected reconstructed events;Raw FT0M status;Events",
               HistType::kTH1D, {AxisSpec{static_cast<int>(FT0MEventStatus::Unavailable), 0.5, static_cast<int>(FT0MEventStatus::Unavailable) + 0.5, "Raw FT0M status"}}, true);
    auto* ft0mStatusAxis = histos.get<TH1>(HIST("hEventFT0MStatus"))->GetXaxis();
    ft0mStatusAxis->SetBinLabel(static_cast<int>(FT0MEventStatus::Valid), "Finite FT0M [0,100)");
    ft0mStatusAxis->SetBinLabel(static_cast<int>(FT0MEventStatus::Negative), "Negative FT0M");
    ft0mStatusAxis->SetBinLabel(static_cast<int>(FT0MEventStatus::OutOfRange), "FT0M >= 100");
    ft0mStatusAxis->SetBinLabel(static_cast<int>(FT0MEventStatus::NonFinite), "Nonfinite FT0M");
    ft0mStatusAxis->SetBinLabel(static_cast<int>(FT0MEventStatus::Unavailable), "FT0M unavailable");

    // Invariant-mass spectra
    histos.add("hInvMassPKUnlike",
               "pK invariant mass (unlike-sign);M_{pK} (GeV/c^{2});Counts",
               HistType::kTH1F, {massAxis});
    histos.add("hInvMassPKLike",
               "pK invariant mass (like-sign);M_{pK} (GeV/c^{2});Counts",
               HistType::kTH1F, {massAxis});

    // THnSparse for invariant-mass analysis (mass, pT, multiplicity/centrality)
    if (lstarEnableSparse.value != 0) {

      if (lstarEnableProxyControls.value != 0) {
        // Axes 3/4/5 partition the SAME inclusive proxy sample in every pair class.
        const std::vector<AxisSpec> controlAxes{
          AxisSpec{400, 1.4, 1.8, "M_{(d/2)K} (GeV/c^{2})"},
          AxisSpec{100, 0., 10., "p_{T}^{(d/2)K} (GeV/c)"}, centAxis,
          AxisSpec{2, -0.5, 1.5, "deuteron hasTOF (0/1)"},
          AxisSpec{2, -0.5, 1.5, "passes full proton PID (0/1)"},
          AxisSpec{2, -0.5, 1.5, "passes proton TPC cut only (0/1)"}};
        histos.add("hLambdaStarProxySelectedTrackControlSparse", "Selected deuteron track controls",
                   HistType::kTHnSparseF,
                   {AxisSpec{100, 0., 10., "p_{T}^{d,full} (GeV/c)"}, controlAxes[3], controlAxes[4], controlAxes[5]}, true);
        auto mapAxes = controlAxes;
        // Full-pK rapidity is recorded, not imposed. Include mass flow bins
        // when projecting these maps back onto the inclusive proxy sample.
        mapAxes.emplace_back(1000, 1.4, 2.4, "M_{pK}^{full} (GeV/c^{2})");
        mapAxes.emplace_back(2, -0.5, 1.5, "passes full-pK rapidity (0/1)");
        for (auto const& pairClass : {"Unlike", "Like", "Mixed"}) {
          if (std::strcmp(pairClass, "Mixed") == 0 && lstarNoMixedEvents.value == 0) {
            continue;
          }
          const std::string spectrumName = std::string("hLambdaStarProxyControl") + pairClass + "Sparse";
          const std::string mapName = std::string("hLambdaStarProxyControlFullPK") + pairClass + "Sparse";
          histos.add(spectrumName.c_str(), "Inclusive proxy acceptance; PID compatibility controls", HistType::kTHnSparseF, controlAxes, true);
          if (lstarEnableFullPKMaps.value != 0) {
            histos.add(mapName.c_str(), "Inclusive proxy acceptance; full-pK control map", HistType::kTHnSparseF, mapAxes, true);
          }
        }
      }

      if (lstarEnableFullPKMaps.value != 0) {
        if (lstarNoMixedEvents.value > 0) {
          histos.add(
            "hLambdaStarProxyVsFullPKMixedSparse",
            "d-selected mixed-event pairs: full-p proton hypothesis vs (d/2)K proxy;"
            "M_{(d/2)K} (GeV/c^{2});"
            "M_{pK}^{full} (GeV/c^{2});"
            "p_{T}^{(d/2)K} (GeV/c);"
            "multiplicity/centrality",
            HistType::kTHnSparseF,
            {AxisSpec{400, 1.4, 1.8,
                      "M_{(d/2)K} (GeV/c^{2})"},

             AxisSpec{400, 1.4, 1.8,
                      "M_{pK}^{full} (GeV/c^{2})"},

             AxisSpec{100, 0., 10.,
                      "p_{T}^{(d/2)K} (GeV/c)"},

             centAxis});
        }
      }

      if (lstarEnableProtonLikeSparse.value != 0) {
        histos.add(
          "hLambdaStarProxyProtonLikeUnlikeSparse",
          "d-selected and proton-like OS proxy candidates;"
          "M_{(d/2)K} (GeV/c^{2});"
          "p_{T}^{(d/2)K} (GeV/c);"
          "multiplicity/centrality",
          HistType::kTHnSparseF,
          {AxisSpec{400, 1.4, 1.8,
                    "M_{(d/2)K} (GeV/c^{2})"},

           AxisSpec{100, 0., 10.,
                    "p_{T}^{(d/2)K} (GeV/c)"},

           centAxis});
      }

      if (lstarEnableProtonLikeSparse.value != 0) {
        histos.add(
          "hLambdaStarProxyProtonLikeLikeSparse",
          "d-selected and proton-like LS proxy candidates;"
          "M_{(d/2)K} (GeV/c^{2});"
          "p_{T}^{(d/2)K} (GeV/c);"
          "multiplicity/centrality",
          HistType::kTHnSparseF,
          {AxisSpec{400, 1.4, 1.8,
                    "M_{(d/2)K} (GeV/c^{2})"},

           AxisSpec{100, 0., 10.,
                    "p_{T}^{(d/2)K} (GeV/c)"},

           centAxis});
      }

      histos.add("hLambdaStarPKUnlikeSparse",
                 "#Lambda^{*}(1520) pK unlike-sign candidates;M_{pK} (GeV/c^{2});p_{T}^{pK} (GeV/c);multiplicity/centrality",
                 HistType::kTHnSparseF,
                 {AxisSpec{400, 1.4, 1.8, "M_{pK} (GeV/c^{2})"},
                  AxisSpec{100, 0., 10., "p_{T}^{pK} (GeV/c)"}, centAxis});

      histos.add("hLambdaStarPKLikeSparse",
                 "#Lambda^{*}(1520) pK like-sign candidates;M_{pK} (GeV/c^{2});p_{T}^{pK} (GeV/c);multiplicity/centrality",
                 HistType::kTHnSparseF,
                 {AxisSpec{400, 1.4, 1.8, "M_{pK} (GeV/c^{2})"},
                  AxisSpec{100, 0., 10., "p_{T}^{pK} (GeV/c)"}, centAxis});

      if (lstarNoMixedEvents.value > 0) {
        histos.add("hLambdaStarPKMixedSparse",
                   "#Lambda^{*}(1520) pK mixed-event candidates;M_{pK} (GeV/c^{2});p_{T}^{pK} (GeV/c);multiplicity/centrality",
                   HistType::kTHnSparseF,
                   {AxisSpec{400, 1.4, 1.8, "M_{pK} (GeV/c^{2})"},
                    AxisSpec{100, 0., 10., "p_{T}^{pK} (GeV/c)"}, centAxis});
      }

      // THnSparse for deuteron-proxy invariant-mass analysis (mass, pT, multiplicity/centrality)
      histos.add("hLambdaStarProxySparse",
                 "#Lambda^{*}(1520) deuteron-proxy candidates;M_{p_{proxy}K} (GeV/c^{2});p_{T}^{p_{proxy}K} (GeV/c);multiplicity/centrality",
                 HistType::kTHnSparseF,
                 {AxisSpec{400, 1.4, 1.8, "M_{p_{proxy}K} (GeV/c^{2})"},
                  AxisSpec{100, 0., 10., "p_{T}^{p_{proxy}K} (GeV/c)"}, centAxis});

      histos.add("hLambdaStarProxyLikeSparse",
                 "#Lambda^{*}(1520) deuteron-proxy like-sign candidates;M_{p_{proxy}K} (GeV/c^{2});p_{T}^{p_{proxy}K} (GeV/c);multiplicity/centrality",
                 HistType::kTHnSparseF,
                 {AxisSpec{400, 1.4, 1.8, "M_{p_{proxy}K} (GeV/c^{2})"},
                  AxisSpec{100, 0., 10., "p_{T}^{p_{proxy}K} (GeV/c)"}, centAxis});

      if (lstarNoMixedEvents.value > 0) {
        histos.add("hLambdaStarProxyMixedSparse",
                   "#Lambda^{*}(1520) deuteron-proxy mixed-event candidates;M_{p_{proxy}K} (GeV/c^{2});p_{T}^{p_{proxy}K} (GeV/c);multiplicity/centrality",
                   HistType::kTHnSparseF,
                   {AxisSpec{400, 1.4, 1.8, "M_{p_{proxy}K} (GeV/c^{2})"},
                    AxisSpec{100, 0., 10., "p_{T}^{p_{proxy}K} (GeV/c)"}, centAxis});
      }

      if (lstarEnableFullPKMaps.value != 0) {
        histos.add(
          "hLambdaStarProxyVsFullPKUnlikeSparse",
          "d-selected OS pairs: full-p proton hypothesis vs (d/2)K proxy;"
          "M_{(d/2)K} (GeV/c^{2});"
          "M_{pK}^{full} (GeV/c^{2});"
          "p_{T}^{(d/2)K} (GeV/c);"
          "multiplicity/centrality",
          HistType::kTHnSparseF,
          {AxisSpec{400, 1.4, 1.8,
                    "M_{(d/2)K} (GeV/c^{2})"},

           AxisSpec{400, 1.4, 1.8,
                    "M_{pK}^{full} (GeV/c^{2})"},

           AxisSpec{100, 0., 10.,
                    "p_{T}^{(d/2)K} (GeV/c)"},

           centAxis});
      }

      if (lstarEnableFullPKMaps.value != 0) {
        histos.add(
          "hLambdaStarProxyVsFullPKLikeSparse",
          "d-selected LS pairs: full-p proton hypothesis vs (d/2)K proxy;"
          "M_{(d/2)K} (GeV/c^{2});"
          "M_{pK}^{full} (GeV/c^{2});"
          "p_{T}^{(d/2)K} (GeV/c);"
          "multiplicity/centrality",
          HistType::kTHnSparseF,
          {AxisSpec{400, 1.4, 1.8,
                    "M_{(d/2)K} (GeV/c^{2})"},

           AxisSpec{400, 1.4, 1.8,
                    "M_{pK}^{full} (GeV/c^{2})"},

           AxisSpec{100, 0., 10.,
                    "p_{T}^{(d/2)K} (GeV/c)"},

           centAxis});
      }
    }

    // Deuteron-proxy invariant mass (p_{proxy} from d/2 combined with K)
    histos.add("hDeuteronProxyMass",
               "#Lambda^{*} proxy invariant mass from (d/2 + K);M_{p_{proxy}K} (GeV/c^{2});Counts",
               HistType::kTH1F, {massAxis});

    if (lstarEnablePidQA.value != 0) {
      // Inclusive PID QA before final candidate PID cuts: after track-quality cuts only
      // TPC dE/dx vs total momentum
      histos.add("hTPCdEdxVsP",
                 "TPC dE/dx vs p;p (GeV/c);dE/dx (arb. units);Counts",
                 HistType::kTH2F, {pAxis, dEdxAxis});

      // TOF #beta vs total momentum
      histos.add("hTOFBetaVsP",
                 "TOF #beta vs p;p (GeV/c);#beta_{TOF};Counts",
                 HistType::kTH2F, {pAxis, betaAxis});

      histos.add("hHasTOFVsP",
                 "TOF matching flag vs p;p (GeV/c);has TOF match;Counts",
                 HistType::kTH2F, {pAxis, tofMatchAxis});

      histos.add("hHasTOFVsPt",
                 "TOF matching flag vs p_{T};p_{T} (GeV/c);has TOF match;Counts",
                 HistType::kTH2F, {ptAxis, tofMatchAxis});

      // --- Per-species inclusive PID QA before final candidate PID cuts ---
      // Species tagging here uses classifyPidSpecies() and is meant only for QA.
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

      // Inclusive tagged nSigma QA before final candidate PID cuts: after track-quality cuts only.
      // These are separate from the final-candidate PID QA histograms.
      histos.add("hNsigmaTPCPionTaggedVsP",
                 "TPC n#sigma_{#pi} vs p for tagged #pi;p (GeV/c);n#sigma^{TPC}_{#pi};Counts",
                 HistType::kTH2F, {pAxis, nsAxis});
      histos.add("hNsigmaTOFPionTaggedVsP",
                 "TOF n#sigma_{#pi} vs p for tagged #pi;p (GeV/c);n#sigma^{TOF}_{#pi};Counts",
                 HistType::kTH2F, {pAxis, nsAxis});

      histos.add("hNsigmaTPCKaonTaggedVsP",
                 "TPC n#sigma_{K} vs p for tagged K;p (GeV/c);n#sigma^{TPC}_{K};Counts",
                 HistType::kTH2F, {pAxis, nsAxis});
      histos.add("hNsigmaTOFKaonTaggedVsP",
                 "TOF n#sigma_{K} vs p for tagged K;p (GeV/c);n#sigma^{TOF}_{K};Counts",
                 HistType::kTH2F, {pAxis, nsAxis});

      histos.add("hNsigmaTPCProtonTaggedVsP",
                 "TPC n#sigma_{p} vs p for tagged p;p (GeV/c);n#sigma^{TPC}_{p};Counts",
                 HistType::kTH2F, {pAxis, nsAxis});
      histos.add("hNsigmaTOFProtonTaggedVsP",
                 "TOF n#sigma_{p} vs p for tagged p;p (GeV/c);n#sigma^{TOF}_{p};Counts",
                 HistType::kTH2F, {pAxis, nsAxis});

      histos.add("hNsigmaTPCDeuteronTaggedVsP",
                 "TPC n#sigma_{d} vs p for tagged d;p (GeV/c);n#sigma^{TPC}_{d};Counts",
                 HistType::kTH2F, {pAxis, nsAxis});
      histos.add("hNsigmaTOFDeuteronTaggedVsP",
                 "TOF n#sigma_{d} vs p for tagged d;p (GeV/c);n#sigma^{TOF}_{d};Counts",
                 HistType::kTH2F, {pAxis, nsAxis});
    }

    if (doprocessMCQA && lstarFillMCTruth.value != 0) {
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
    }

    if (lstarEnablePidQA.value != 0) {
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

      // --- Candidate PID QA after final selected K/p/d PID cuts ---
      // These histograms are needed for PID studies in the same pT intervals used by the analysis.
      histos.add("hTOFBetaVsPt_K",
                 "TOF #beta vs p_{T} for selected K;p_{T} (GeV/c);#beta_{TOF};Counts",
                 HistType::kTH2F, {ptAxis, betaAxis});
      histos.add("hTOFBetaVsPt_P",
                 "TOF #beta vs p_{T} for selected p;p_{T} (GeV/c);#beta_{TOF};Counts",
                 HistType::kTH2F, {ptAxis, betaAxis});
      histos.add("hTOFBetaVsPt_D",
                 "TOF #beta vs p_{T} for selected d;p_{T} (GeV/c);#beta_{TOF};Counts",
                 HistType::kTH2F, {ptAxis, betaAxis});

      histos.add("hNsigmaTPCKaonVsPt",
                 "TPC n#sigma_{K} vs p_{T};p_{T} (GeV/c);n#sigma^{TPC}_{K};Counts",
                 HistType::kTH2F, {ptAxis, nsAxis});
      histos.add("hNsigmaTOFKaonVsPt",
                 "TOF n#sigma_{K} vs p_{T};p_{T} (GeV/c);n#sigma^{TOF}_{K};Counts",
                 HistType::kTH2F, {ptAxis, nsAxis});

      histos.add("hNsigmaTPCProtonVsP",
                 "TPC n#sigma_{p} vs p;p (GeV/c);n#sigma^{TPC}_{p};Counts",
                 HistType::kTH2F, {pAxis, nsAxis});
      histos.add("hNsigmaTOFProtonVsP",
                 "TOF n#sigma_{p} vs p;p (GeV/c);n#sigma^{TOF}_{p};Counts",
                 HistType::kTH2F, {pAxis, nsAxis});
      histos.add("hNsigmaTPCProtonVsPt",
                 "TPC n#sigma_{p} vs p_{T};p_{T} (GeV/c);n#sigma^{TPC}_{p};Counts",
                 HistType::kTH2F, {ptAxis, nsAxis});
      histos.add("hNsigmaTOFProtonVsPt",
                 "TOF n#sigma_{p} vs p_{T};p_{T} (GeV/c);n#sigma^{TOF}_{p};Counts",
                 HistType::kTH2F, {ptAxis, nsAxis});

      histos.add("hNsigmaTPCDeuteronVsPt",
                 "TPC n#sigma_{d} vs p_{T};p_{T} (GeV/c);n#sigma^{TPC}_{d};Counts",
                 HistType::kTH2F, {ptAxis, nsAxis});
      histos.add("hNsigmaTOFDeuteronVsPt",
                 "TOF n#sigma_{d} vs p_{T};p_{T} (GeV/c);n#sigma^{TOF}_{d};Counts",
                 HistType::kTH2F, {ptAxis, nsAxis});
      histos.add(
        "hNsigmaTPCProtonForSelectedDeuteronVsPt",
        "TPC n#sigma_{p} for selected d candidates vs p_{T};"
        "p_{T} (GeV/c);"
        "n#sigma^{TPC}_{p};Counts",
        HistType::kTH2F,
        {ptAxis, nsAxis});

      histos.add(
        "hNsigmaTOFProtonForSelectedDeuteronVsPt",
        "TOF n#sigma_{p} for selected d candidates vs p_{T};"
        "p_{T} (GeV/c);"
        "n#sigma^{TOF}_{p};Counts",
        HistType::kTH2F,
        {ptAxis, nsAxis});

      histos.add(
        "hTPCNsigmaDeVsPrForSelectedDeuteron",
        "TPC PID overlap for selected d candidates;"
        "n#sigma^{TPC}_{d};"
        "n#sigma^{TPC}_{p};Counts",
        HistType::kTH2F,
        {nsAxis, nsAxis});

      histos.add(
        "hTOFNsigmaDeVsPrForSelectedDeuteron",
        "TOF PID overlap for selected d candidates;"
        "n#sigma^{TOF}_{d};"
        "n#sigma^{TOF}_{p};Counts",
        HistType::kTH2F,
        {nsAxis, nsAxis});

      histos.add("hTPCvsTOFNsigma_K",
                 "TPC vs TOF n#sigma for selected K;n#sigma^{TPC}_{K};n#sigma^{TOF}_{K};Counts",
                 HistType::kTH2F, {nsAxis, nsAxis});
      histos.add("hTPCvsTOFNsigma_P",
                 "TPC vs TOF n#sigma for selected p;n#sigma^{TPC}_{p};n#sigma^{TOF}_{p};Counts",
                 HistType::kTH2F, {nsAxis, nsAxis});
      histos.add("hTPCvsTOFNsigma_D",
                 "TPC vs TOF n#sigma for selected d;n#sigma^{TPC}_{d};n#sigma^{TOF}_{d};Counts",
                 HistType::kTH2F, {nsAxis, nsAxis});

      // --- DCA QA for final selected K/p/d candidates ---
      // Filled only when lstarEnablePidQA = 1.
      histos.add("hDCAxyVsPt_K",
                 "DCA_{xy} vs p_{T} distribution for selected K;p_{T} (GeV/c);DCA_{xy} (cm);Counts",
                 HistType::kTH2F, {ptAxis, dcaXYAxis});
      histos.add("hDCAzVsPt_K",
                 "DCA_{z} vs p_{T} distribution for selected K;p_{T} (GeV/c);DCA_{z} (cm);Counts",
                 HistType::kTH2F, {ptAxis, dcaZAxis});

      histos.add("hDCAxyVsPt_P",
                 "DCA_{xy} vs p_{T} distribution for selected p;p_{T} (GeV/c);DCA_{xy} (cm);Counts",
                 HistType::kTH2F, {ptAxis, dcaXYAxis});
      histos.add("hDCAzVsPt_P",
                 "DCA_{z} vs p_{T} distribution for selected p;p_{T} (GeV/c);DCA_{z} (cm);Counts",
                 HistType::kTH2F, {ptAxis, dcaZAxis});

      histos.add("hDCAxyVsPt_D",
                 "DCA_{xy} vs p_{T} distribution for selected d;p_{T} (GeV/c);DCA_{xy} (cm);Counts",
                 HistType::kTH2F, {ptAxis, dcaXYAxis});
      histos.add("hDCAzVsPt_D",
                 "DCA_{z} vs p_{T} distribution for selected d;p_{T} (GeV/c);DCA_{z} (cm);Counts",
                 HistType::kTH2F, {ptAxis, dcaZAxis});
    }
    bookMC();
    bookProxyDiagnostics();
  }

  static constexpr int PdgLambdaStar = static_cast<int>(PDG_t::kLambda1520);
  static constexpr int PdgDeuteron = static_cast<int>(o2::constants::physics::Pdg::kDeuteron);
  // MC outputs are 2D mass/pT histograms. They add no high-dimensional maps.
  std::array<std::shared_ptr<TH2>, 5> mcGenSpectra{};
  std::array<std::array<std::shared_ptr<TH2>, 2>, 4> mcRecoProxy{};
  std::shared_ptr<TH2> mcRecoPKTrue, mcPKResponse, mcProxyResponse;

  // Fixed bounded diagnostics: the four analysis pT bins plus ROOT flow bins.
  // Existing native 1-MeV/0.1-GeV histograms are not changed or replaced.
  static constexpr std::array<double, 5> DiagnosticPairPtEdges{0.6, 1.3, 1.4, 1.9, 3.6};
  static constexpr std::array<double, 9> ProbeMomentumEdges{0., 0.5, 0.8, 1.2, 1.8, 2.5, 3.5, 5., 10.};
  static constexpr std::array<const char*, 7> SpeciesNames{"e", "mu", "pi", "K", "p", "d", "remaining"};
  static constexpr int SpeciesCount = SpeciesNames.size();
  static constexpr int SpeciesPairs = SpeciesCount * SpeciesCount;
  static constexpr int UnmatchedSpeciesTag = SpeciesPairs;
  static constexpr int ChargeMismatchSpeciesTag = SpeciesPairs + 1;
  static constexpr int SpeciesTagCount = SpeciesPairs + 2;
  static constexpr int ProtonSpecies = 4;
  static constexpr int DeuteronSpecies = 5;
  static constexpr int KaonSpecies = 3;
  static constexpr int SelectedPIDStates = 16;
  static constexpr int ProbeRegions = 4;
  static constexpr int ProbePairPtBins = DiagnosticPairPtEdges.size() - 1;
  static constexpr int ProbeMomentumBins = ProbeMomentumEdges.size(); // includes p >= 10
  static constexpr int ProbeTagCount = ProbePairPtBins * ProbeMomentumBins * ProbeRegions;
  static constexpr int ParentClasses = 10;
  static constexpr int DiagnosticMassBins = 200;
  static constexpr int ProbeMassBins = 100;
  static constexpr float ProbeTPCMax = 6.f;
  static constexpr float ProbeTOFMax = 8.f;
  static constexpr double MultMin = 0.;
  static constexpr double MultMax = 100.;
  static constexpr std::array<const char*, 2> PairSigns{"Unlike", "Like"};
  static constexpr std::array<const char*, 2> ProbeRoles{"De", "Ka"};
  static constexpr std::array<const char*, ParentClasses> ParentNames{
    "no retained common mother", "K0S", "Lambda1115", "Kstar", "phi", "Lambda1520",
    "other retained common mother", "multiple retained common mothers", "same MC particle", "invalid truth labels"};

  std::array<std::shared_ptr<TH3>, 2> mcProxySpecies{};
  std::array<std::shared_ptr<TH3>, 2> mcProxyOtherParents{};
  std::array<std::array<std::shared_ptr<TH3>, 2>, 2> probeMass2{}, probeTPC{}, probeTOF{};
  std::array<std::shared_ptr<TH3>, 2> probePairRegions{}, mcProbeSpecies{};
  std::shared_ptr<TH2> probeTrackStatus;
  std::shared_ptr<TH1> probeEventScope;

  struct PIDProbeCand {
    float px = 0.f, py = 0.f, pz = 0.f;
    double momentum = 0., mass2 = 0.;
    float nsTPC = 0.f, nsTOF = 0.f;
    bool nominal = false;
    int charge = 0;
    int64_t tid = -1;
    TruthInfo truth{};
  };

  static bool regularMultiplicity(float mult)
  {
    return std::isfinite(mult) && mult >= MultMin && mult < MultMax;
  }

  static int speciesGroup(int pdg)
  {
    switch (std::abs(pdg)) {
      case PdgElectron:
        return 0;
      case PdgMuon:
        return 1;
      case PdgPion:
        return 2;
      case PdgKaon:
        return KaonSpecies;
      case PdgProton:
        return ProtonSpecies;
      case PdgDeuteron:
        return DeuteronSpecies;
      default:
        return SpeciesCount - 1;
    }
  }

  int speciesPairTag(TruthInfo const& a, TruthInfo const& b, int chargeA, int chargeB)
  {
    if (a.id < 0 || b.id < 0) {
      return UnmatchedSpeciesTag;
    }
    if (chargeFromPdg(a.pdg) != chargeA || chargeFromPdg(b.pdg) != chargeB) {
      return ChargeMismatchSpeciesTag;
    }
    return SpeciesCount * speciesGroup(a.pdg) + speciesGroup(b.pdg);
  }

  static bool otherSpeciesTag(int tag)
  {
    return tag != UnmatchedSpeciesTag && tag != SpeciesCount * DeuteronSpecies + KaonSpecies && tag != SpeciesCount * ProtonSpecies + KaonSpecies;
  }

  static int retainedParentClass(TruthInfo const& a, TruthInfo const& b)
  {
    if (a.id < 0 || b.id < 0) {
      return ParentClasses - 1;
    }
    if (a.id == b.id) {
      return ParentClasses - 2;
    }
    int64_t commonId = -1;
    int commonPdg = 0;
    for (auto const& ma : a.diagnosticMothers) {
      for (auto const& mb : b.diagnosticMothers) {
        if (ma.id != mb.id) {
          continue;
        }
        if (commonId >= 0 && commonId != ma.id) {
          return ParentClasses - 3;
        }
        commonId = ma.id;
        commonPdg = std::abs(ma.pdg);
      }
    }
    if (commonId < 0) {
      return 0; // absence of a retained link does not establish absence of ancestry
    }
    if (commonPdg == static_cast<int>(PDG_t::kK0Short)) {
      return 1;
    }
    if (commonPdg == static_cast<int>(PDG_t::kLambda0)) {
      return 2;
    }
    if (commonPdg == static_cast<int>(o2::constants::physics::Pdg::kK0Star892) || commonPdg == static_cast<int>(o2::constants::physics::Pdg::kKPlusStar892)) {
      return 3;
    }
    if (commonPdg == static_cast<int>(o2::constants::physics::Pdg::kPhi)) {
      return 4;
    }
    if (commonPdg == PdgLambdaStar) {
      return 5;
    }
    return 6;
  }

  static int probePairPtBin(float pt)
  {
    if (!std::isfinite(pt) || pt < DiagnosticPairPtEdges.front() || pt >= DiagnosticPairPtEdges.back()) {
      return -1;
    }
    return static_cast<int>(std::upper_bound(DiagnosticPairPtEdges.begin(), DiagnosticPairPtEdges.end(), pt) - DiagnosticPairPtEdges.begin()) - 1;
  }

  static int probeMomentumBin(double momentum)
  {
    return std::clamp(static_cast<int>(std::upper_bound(ProbeMomentumEdges.begin(), ProbeMomentumEdges.end(), momentum) - ProbeMomentumEdges.begin()) - 1, 0, ProbeMomentumBins - 1);
  }

  static int probeTag(int pairBin, double momentum, int region)
  {
    return pairBin + ProbePairPtBins * (probeMomentumBin(momentum) + ProbeMomentumBins * region);
  }

  static std::string speciesTagName(int tag)
  {
    if (tag == UnmatchedSpeciesTag) {
      return "unmatched";
    }
    if (tag == ChargeMismatchSpeciesTag) {
      return "charge-inconsistent";
    }
    return std::string(SpeciesNames[tag / SpeciesCount]) + ":" + SpeciesNames[tag % SpeciesCount];
  }

  void bookProxyDiagnostics()
  {
    if (lstarEnableMCProxySpecies.value != 0 && !doprocessMCReco) {
      LOGF(fatal, "lstarEnableMCProxySpecies requires processMCReco");
    }
    if (lstarEnableProxyPIDProbes.value != 0) {
      if (!doprocess && !doprocessMCReco) {
        LOGF(fatal, "PID probes require the data or reconstructed-MC process");
      }
      for (auto const& cut : {lstarProxyProbeTPCDeMax.value, lstarProxyProbeTPCKaMax.value}) {
        if (!std::isfinite(cut) || cut <= 0.f || cut > ProbeTPCMax) {
          LOGF(fatal, "Probe TPC limits must lie in (0,6]");
        }
      }
      for (auto const& cut : {lstarProxyProbeTOFDeMax.value, lstarProxyProbeTOFKaMax.value}) {
        if (!std::isfinite(cut) || cut <= 0.f || cut > ProbeTOFMax) {
          LOGF(fatal, "Probe TOF limits must lie in (0,8]");
        }
      }
    }
    const AxisSpec pairPt{std::vector<double>(DiagnosticPairPtEdges.begin(), DiagnosticPairPtEdges.end()), "Proxy pair pT (GeV/c)"};
    const AxisSpec mcMass{DiagnosticMassBins, 1.4, 1.8, "M_(d/2)K (GeV/c2)"};
    const AxisSpec probeMass{ProbeMassBins, 1.4, 1.8, "M_(d/2)K (GeV/c2)"};
    uint64_t payload = 0;
    auto add3 = [&](std::string const& name, std::string const& title, std::vector<AxisSpec> const& axes) {
      auto h = histos.add<TH3>(name.c_str(), title.c_str(), HistType::kTH3D, axes, true);
      payload += static_cast<uint64_t>(h->GetNcells()) * 2 * sizeof(double);
      return h;
    };
    for (size_t sign = 0; sign < PairSigns.size(); ++sign) {
      const std::string suffix = PairSigns[sign];
      if (lstarEnableMCProxySpecies.value != 0) {
        mcProxySpecies[sign] = add3("hMCProxySpecies" + suffix, "Selected proxy;species+51*(dTOF+2*kTOF+4*pTPC+8*regularMult)",
                                    {mcMass, pairPt, AxisSpec{SpeciesTagCount * SelectedPIDStates, -0.5, SpeciesTagCount * SelectedPIDStates - 0.5, "Species and detector state"}});
        for (int state = 0; state < SelectedPIDStates; ++state) {
          for (int tag = 0; tag < SpeciesTagCount; ++tag) {
            const std::string label = speciesTagName(tag) + ";state=" + std::to_string(state);
            mcProxySpecies[sign]->GetZaxis()->SetBinLabel(tag + SpeciesTagCount * state + 1, label.c_str());
          }
        }
        mcProxyOtherParents[sign] = add3("hMCProxyOtherParents" + suffix, "Selected Other: retained immediate common-mother class only",
                                         {mcMass, pairPt, AxisSpec{ParentClasses, -0.5, ParentClasses - 0.5, "Retained parent class"}});
        for (int parent = 0; parent < ParentClasses; ++parent) {
          mcProxyOtherParents[sign]->GetZaxis()->SetBinLabel(parent + 1, ParentNames[parent]);
        }
      }
      if (lstarEnableProxyPIDProbes.value == 0) {
        continue;
      }
      probePairRegions[sign] = add3("hProxyPIDProbeRegions" + suffix, "Loose matched PID probes;regular multiplicity only;region=dNominal+2*kNominal",
                                    {probeMass, pairPt, AxisSpec{ProbeRegions, -0.5, ProbeRegions - 0.5, "Nominal-selection region"}});
      for (size_t role = 0; role < ProbeRoles.size(); ++role) {
        const std::string stem = std::string("hProxyPIDProbe") + ProbeRoles[role];
        const AxisSpec m2 = role == 0 ? AxisSpec{120, -1., 7., "Full-track TOF mass2 (GeV2/c4)"} : AxisSpec{80, -1., 2., "Full-track TOF mass2 (GeV2/c4)"};
        const AxisSpec packed{ProbeTagCount, -0.5, ProbeTagCount - 0.5, "pairBin+4*(momentumBin+9*region)"};
        probeMass2[role][sign] = add3(std::string(stem).append("Mass2").append(suffix), "Loose both-TOF proxy probes;TPC/TOF envelopes;regular multiplicity", {probeMass, m2, packed});
        probeTPC[role][sign] = add3(std::string(stem).append("TPC").append(suffix), "Same loose both-TOF proxy probes;TPC nSigma;regular multiplicity", {probeMass, AxisSpec{60, -ProbeTPCMax, ProbeTPCMax, "TPC nSigma for slot hypothesis"}, packed});
        probeTOF[role][sign] = add3(std::string(stem).append("TOF").append(suffix), "Same loose both-TOF proxy probes;TOF nSigma;regular multiplicity", {probeMass, AxisSpec{80, -ProbeTOFMax, ProbeTOFMax, "TOF nSigma for slot hypothesis"}, packed});
        for (int region = 0; region < ProbeRegions; ++region) {
          for (int momentum = 0; momentum < ProbeMomentumBins; ++momentum) {
            for (int pair = 0; pair < ProbePairPtBins; ++pair) {
              const int tag = pair + ProbePairPtBins * (momentum + ProbeMomentumBins * region);
              const std::string label = "pair=" + std::to_string(pair) + ";p=" + std::to_string(momentum) + ";region=" + std::to_string(region);
              probeMass2[role][sign]->GetZaxis()->SetBinLabel(tag + 1, label.c_str());
              probeTPC[role][sign]->GetZaxis()->SetBinLabel(tag + 1, label.c_str());
              probeTOF[role][sign]->GetZaxis()->SetBinLabel(tag + 1, label.c_str());
            }
          }
        }
      }
      if (doprocessMCReco) {
        mcProbeSpecies[sign] = add3("hMCProxyPIDProbeSpecies" + suffix, "Same loose probe pairs;species+51*nominalRegion;regular multiplicity",
                                    {probeMass, pairPt, AxisSpec{SpeciesTagCount * ProbeRegions, -0.5, SpeciesTagCount * ProbeRegions - 0.5, "Species plus nominal region"}});
        for (int region = 0; region < ProbeRegions; ++region) {
          for (int tag = 0; tag < SpeciesTagCount; ++tag) {
            const std::string label = speciesTagName(tag) + ";region=" + std::to_string(region);
            mcProbeSpecies[sign]->GetZaxis()->SetBinLabel(tag + SpeciesTagCount * region + 1, label.c_str());
          }
        }
      }
    }
    if (lstarEnableProxyPIDProbes.value != 0) {
      probeTrackStatus = histos.add<TH2>("hProxyPIDProbeTrackStatus", "After quality and loose TPC;slot;stage", HistType::kTH2D,
                                         {AxisSpec{2, -0.5, 1.5, "Slot: d=0,K=1"}, AxisSpec{5, -0.5, 4.5, "Stage"}}, true);
      const std::array<const char*, 5> labels{"loose TPC", "no TOF", "invalid TOF response", "outside loose TOF", "accepted probe track"};
      for (size_t i = 0; i < labels.size(); ++i) {
        probeTrackStatus->GetYaxis()->SetBinLabel(i + 1, labels[i]);
      }
      probeEventScope = histos.add<TH1>("hProxyPIDProbeEventScope", "Selected events;regular multiplicity flag;events", HistType::kTH1D, {AxisSpec{2, -0.5, 1.5, "0=outside,1=regular [0,100)"}}, true);
    }
    if (payload > 0) {
      LOGF(info, "Proxy diagnostic TH3 content+Sumw2 payload: %.2f MiB per wagon (fixed bins; object/axis overhead extra)", static_cast<double>(payload) / (1024. * 1024.));
    }
  }

  void fillMCProxySpecies(ProxyCand const& d, KaonCand const& k, double mass, float pt, float mult)
  {
    const int tag = speciesPairTag(d.truth, k.truth, d.charge, k.charge);
    const int state = static_cast<int>(d.hasTOF) + 2 * static_cast<int>(k.hasTOF) + 4 * static_cast<int>(d.protonTPCCompatible) + 8 * static_cast<int>(regularMultiplicity(mult));
    const size_t sign = d.charge * k.charge < 0 ? 0 : 1;
    mcProxySpecies[sign]->Fill(mass, pt, tag + SpeciesTagCount * state);
    if (otherSpeciesTag(tag)) {
      mcProxyOtherParents[sign]->Fill(mass, pt, retainedParentClass(d.truth, k.truth));
    }
  }

  template <bool IsMC, typename TTracks>
  void fillProxyPIDProbes(TTracks const& tracks, float mult, int collisionId)
  {
    // This separate path must not feed nominal candidate vectors or mixing pools.
    const bool regular = regularMultiplicity(mult);
    probeEventScope->Fill(regular);
    if (!regular) {
      return;
    }
    std::array<std::vector<PIDProbeCand>, 2> candidates;
    for (auto const& track : tracks) {
      if (track.pt() < lstarCutPtMin.value || std::abs(track.eta()) > lstarCutEtaMax.value || !passTrackQuality(track) || track.sign() == 0) {
        continue;
      }
      const std::array<float, 2> tpc{track.tpcNSigmaDe(), track.tpcNSigmaKa()};
      const std::array<float, 2> tof{track.tofNSigmaDe(), track.tofNSigmaKa()};
      const std::array<float, 2> tpcMax{lstarProxyProbeTPCDeMax.value, lstarProxyProbeTPCKaMax.value};
      const std::array<float, 2> tofMax{lstarProxyProbeTOFDeMax.value, lstarProxyProbeTOFKaMax.value};
      for (size_t role = 0; role < ProbeRoles.size(); ++role) {
        if (!std::isfinite(tpc[role]) || std::abs(tpc[role]) >= tpcMax[role]) {
          continue;
        }
        probeTrackStatus->Fill(role, 0);
        if (!track.hasTOF()) {
          probeTrackStatus->Fill(role, 1);
          continue;
        }
        const double beta = track.beta();
        const double momentum = static_cast<double>(track.pt()) * std::cosh(static_cast<double>(track.eta()));
        if (!std::isfinite(beta) || beta <= TofBetaMin || beta >= TofBetaMax || !std::isfinite(momentum) || momentum <= 0. || !std::isfinite(tof[role])) {
          probeTrackStatus->Fill(role, 2);
          continue;
        }
        if (std::abs(tof[role]) >= tofMax[role]) {
          probeTrackStatus->Fill(role, 3);
          continue;
        }
        probeTrackStatus->Fill(role, 4);
        PIDProbeCand probe;
        probe.px = track.pt() * std::cos(track.phi());
        probe.py = track.pt() * std::sin(track.phi());
        probe.pz = track.pt() * std::sinh(track.eta());
        probe.momentum = momentum;
        probe.mass2 = momentum * momentum * (1. / (beta * beta) - 1.);
        probe.nsTPC = tpc[role];
        probe.nsTOF = tof[role];
        probe.charge = static_cast<int>(track.sign());
        probe.tid = track.globalIndex();
        if (role == 0) {
          const bool passProton = passFinalCandidatePID(track.pt(), track.tpcNSigmaPr(), track.tofNSigmaPr(), true,
                                                        lstarCutNsigmaTPCPr.value, lstarCutNsigmaTOFPr.value, lstarPidCircularCutPr.value, lstarPidPtRefPr.value);
          probe.nominal = passFinalCandidatePID(track.pt(), tpc[role], tof[role], true, lstarCutNsigmaTPCDe.value, lstarCutNsigmaTOFDe.value, lstarPidCircularCutDe.value, lstarPidPtRefDe.value, true) &&
                          passOptionalDeuteronTOFExtras(track) && (lstarProxyUseProtonPIDAsDeuteron.value == 0 || passProton);
        } else {
          probe.nominal = passFinalCandidatePID(track.pt(), tpc[role], tof[role], true, lstarCutNsigmaTPCKaon.value, lstarCutNsigmaTOFKaon.value, lstarPidCircularCutKaon.value, lstarPidPtRefKaon.value);
        }
        if constexpr (IsMC) {
          probe.truth = trackTruth(track, collisionId);
        }
        candidates[role].push_back(std::move(probe));
      }
    }
    for (auto const& d : candidates[0]) {
      for (auto const& k : candidates[1]) {
        if (d.tid == k.tid) {
          continue;
        }
        const double mass = invariantMass(ProxyMomentumScale * d.px, ProxyMomentumScale * d.py, ProxyMomentumScale * d.pz, o2::constants::physics::MassProton,
                                          k.px, k.py, k.pz, o2::constants::physics::MassKaonCharged);
        const float px = ProxyMomentumScale * d.px + k.px, py = ProxyMomentumScale * d.py + k.py, pz = ProxyMomentumScale * d.pz + k.pz;
        const float pt = ptFromPxPy(px, py);
        const int pairBin = probePairPtBin(pt);
        const float y = rapidityFromMomentumAndMass(px, py, pz, mass);
        if (pairBin < 0 || !std::isfinite(mass) || !std::isfinite(y) || std::abs(y) > lstarLambdaAbsYMax.value) {
          continue;
        }
        const size_t sign = d.charge * k.charge < 0 ? 0 : 1;
        const int region = static_cast<int>(d.nominal) + 2 * static_cast<int>(k.nominal);
        probePairRegions[sign]->Fill(mass, pt, region);
        probeMass2[0][sign]->Fill(mass, d.mass2, probeTag(pairBin, d.momentum, region));
        probeMass2[1][sign]->Fill(mass, k.mass2, probeTag(pairBin, k.momentum, region));
        probeTPC[0][sign]->Fill(mass, d.nsTPC, probeTag(pairBin, d.momentum, region));
        probeTPC[1][sign]->Fill(mass, k.nsTPC, probeTag(pairBin, k.momentum, region));
        probeTOF[0][sign]->Fill(mass, d.nsTOF, probeTag(pairBin, d.momentum, region));
        probeTOF[1][sign]->Fill(mass, k.nsTOF, probeTag(pairBin, k.momentum, region));
        if constexpr (IsMC) {
          mcProbeSpecies[sign]->Fill(mass, pt, speciesPairTag(d.truth, k.truth, d.charge, k.charge) + SpeciesTagCount * region);
        }
      }
    }
  }

  bool mixingEnabled() const
  {
    return lstarEnableSparse.value != 0 && lstarNoMixedEvents.value > 0;
  }

  void bookMC()
  {
    const std::vector<AxisSpec> massPt{
      AxisSpec{400, 1.4, 1.8, "Invariant mass (GeV/c^{2})"},
      AxisSpec{100, 0., 10., "Pair p_{T} (GeV/c)"}};
    const std::vector<AxisSpec> response{
      AxisSpec{100, 0., 10., "Generated pair p_{T} (GeV/c)"},
      AxisSpec{100, 0., 10., "Reconstructed pair p_{T} (GeV/c)"}};
    if (doprocessMCGen) {
      histos.add("hMCGenEventSelection", "Generated events;Stage;Events", HistType::kTH1D,
                 {AxisSpec{3, 0.5, 3.5, "1=all, 2=vertex, 3=vertex+truth INEL>0"}});
      const std::array<const char*, 5> names{
        "hMCGenPKUnlike", "hMCGenPKLike", "hMCGenProxyUnlike", "hMCGenProxyLike", "hMCGenPKTrueLambdaStar"};
      for (size_t i = 0; i < names.size(); ++i) {
        mcGenSpectra[i] = histos.add<TH2>(names[i], "Generated fiducial pairs", HistType::kTH2D, massPt, true);
      }
    }
    if (doprocessMCReco) {
      const std::array<const char*, 4> categories{"TrueD", "ProtonMisID", "Other", "Unmatched"};
      for (size_t i = 0; i < categories.size(); ++i) {
        for (size_t sign = 0; sign < mcRecoProxy[i].size(); ++sign) {
          const std::string name = std::string("hMCRecoProxy") + categories[i] + (sign == 0 ? "Unlike" : "Like");
          mcRecoProxy[i][sign] = histos.add<TH2>(name.c_str(), "Reconstructed proxy truth component", HistType::kTH2D, massPt, true);
        }
      }
      mcRecoPKTrue = histos.add<TH2>("hMCRecoPKTrueLambdaStar", "Selected true common-parent pK", HistType::kTH2D, massPt, true);
      mcPKResponse = histos.add<TH2>("hMCRecoPKPtResponse", "Selected true common-parent pK;Generated pair pT;Reconstructed pair pT", HistType::kTH2D, response, true);
      mcProxyResponse = histos.add<TH2>("hMCRecoProxyPtResponse", "Selected true primary dK pairs;Generated proxy pair pT;Reconstructed proxy pair pT", HistType::kTH2D, response, true);
    }
  }

  template <typename TParticle>
  static TruthInfo particleTruth(TParticle const& particle, bool retainDiagnosticMothers = false)
  {
    TruthInfo info;
    info.id = particle.globalIndex();
    info.collisionId = particle.mcCollisionId();
    info.pdg = particle.pdgCode();
    info.px = particle.px();
    info.py = particle.py();
    info.pz = particle.pz();
    info.primary = particle.isPhysicalPrimary();
    for (auto const& mother : particle.template mothers_as<aod::McParticles>()) {
      if (retainDiagnosticMothers && mother.mcCollisionId() == info.collisionId) {
        info.diagnosticMothers.push_back({mother.globalIndex(), mother.pdgCode()});
      }
      if (std::abs(mother.pdgCode()) == PdgLambdaStar && mother.mcCollisionId() == info.collisionId) {
        info.lambdaMothers.push_back(mother.globalIndex());
      }
    }
    return info;
  }

  template <typename TTrack>
  static TruthInfo trackTruth(TTrack const& track, int collisionId, bool retainDiagnosticMothers = false)
  {
    if (collisionId < 0 || !track.has_mcParticle()) {
      return {};
    }
    auto particle = track.template mcParticle_as<aod::McParticles>();
    if (particle.mcCollisionId() != collisionId) {
      return {}; // pileup or wrong-collision labels do not enter truth components
    }
    return particleTruth(particle, retainDiagnosticMothers);
  }

  static bool commonLambdaStar(TruthInfo const& proton, TruthInfo const& kaon)
  {
    if (std::abs(proton.pdg) != PdgProton || std::abs(kaon.pdg) != PdgKaon ||
        proton.pdg * kaon.pdg >= 0 || proton.id == kaon.id ||
        proton.collisionId != kaon.collisionId) {
      return false;
    }
    for (auto const& id : proton.lambdaMothers) {
      if (std::find(kaon.lambdaMothers.begin(), kaon.lambdaMothers.end(), id) != kaon.lambdaMothers.end()) {
        return true;
      }
    }
    return false;
  }

  template <typename TCandidate>
  void fillRecoMC(TCandidate const& candidate, KaonCand const& kaon, double mass, float pt, bool proxy)
  {
    const auto& a = candidate.truth;
    const auto& b = kaon.truth;
    if (!proxy) {
      if (commonLambdaStar(a, b) && chargeFromPdg(a.pdg) == candidate.charge && chargeFromPdg(b.pdg) == kaon.charge) {
        mcRecoPKTrue->Fill(mass, pt);
        mcPKResponse->Fill(ptFromPxPy(a.px + b.px, a.py + b.py), pt);
      }
      return;
    }
    size_t category = 2; // all other labelled species, secondaries included
    if (a.id < 0 || b.id < 0) {
      category = 3;
    } else if (std::abs(b.pdg) == PdgKaon && chargeFromPdg(b.pdg) == kaon.charge && chargeFromPdg(a.pdg) == candidate.charge) {
      if (std::abs(a.pdg) == PdgDeuteron) {
        category = 0;
      } else if (std::abs(a.pdg) == PdgProton) {
        category = 1;
      }
    }
    const size_t sign = candidate.charge * kaon.charge < 0 ? 0 : 1;
    mcRecoProxy[category][sign]->Fill(mass, pt);
    if (category == 0 && a.primary && b.primary && a.id != b.id) {
      mcProxyResponse->Fill(ptFromPxPy(ProxyMomentumScale * a.px + b.px, ProxyMomentumScale * a.py + b.py), pt);
    }
    // TrueD means true d + true K, not a truth tag of resonance-fed coalescence.
  }

  template <typename TCollision, typename TParticles>
  void fillGenerated(TCollision const& collision, TParticles const& particles)
  {
    histos.fill(HIST("hMCGenEventSelection"), 1.);
    if (!(std::abs(collision.posZ()) < lstarCutVertex.value)) {
      return;
    }
    histos.fill(HIST("hMCGenEventSelection"), 2.);
    // Truth INEL>0 uses charged physical primaries in |eta|<1, with no pT cut.
    // It is independent of detector trigger/vertex reconstruction efficiency.
    bool inel = false;
    for (auto const& particle : particles) {
      if (particle.isPhysicalPrimary() && std::abs(particle.eta()) < 1.f && chargeFromPdg(particle.pdgCode()) != 0) {
        inel = true;
        break;
      }
    }
    if (lstarRequireINELgt0.value != 0 && !inel) {
      return;
    }
    histos.fill(HIST("hMCGenEventSelection"), 3.);
    std::vector<TruthInfo> protons, kaons, deuterons;
    for (auto const& particle : particles) {
      // Use the DATA daughter kinematic acceptance for generated pair spectra.
      // lstarCutMC* remains specific to the optional legacy track-QA process.
      if (!(particle.pt() > lstarCutPtMin.value) || !(std::abs(particle.eta()) < lstarCutEtaMax.value) ||
          (lstarMCPrimaryOnly.value != 0 && !particle.isPhysicalPrimary())) {
        continue;
      }
      switch (std::abs(particle.pdgCode())) {
        case PdgProton:
          protons.push_back(particleTruth(particle));
          break;
        case PdgKaon:
          kaons.push_back(particleTruth(particle));
          break;
        case PdgDeuteron:
          deuterons.push_back(particleTruth(particle));
          break;
        default:
          break;
      }
    }
    auto fillPairs = [&](auto const& candidates, bool proxy) {
      const float scale = proxy ? ProxyMomentumScale : 1.f;
      for (auto const& a : candidates) {
        for (auto const& b : kaons) {
          if (a.id == b.id) {
            continue;
          }
          const double mass = invariantMass(scale * a.px, scale * a.py, scale * a.pz, o2::constants::physics::MassProton,
                                            b.px, b.py, b.pz, o2::constants::physics::MassKaonCharged);
          const float px = scale * a.px + b.px, py = scale * a.py + b.py, pz = scale * a.pz + b.pz;
          if (std::abs(rapidityFromMomentumAndMass(px, py, pz, mass)) > lstarLambdaAbsYMax.value) {
            continue;
          }
          const float pt = ptFromPxPy(px, py);
          const size_t sign = chargeFromPdg(a.pdg) * chargeFromPdg(b.pdg) < 0 ? 0 : 1;
          mcGenSpectra[(proxy ? 2 : 0) + sign]->Fill(mass, pt);
          if (!proxy && commonLambdaStar(a, b)) {
            mcGenSpectra[4]->Fill(mass, pt);
          }
        }
      }
    };
    fillPairs(protons, false);
    fillPairs(deuterons, true);
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

  template <typename TTrack>
  bool passOptionalDeuteronTOFExtras(const TTrack& trk) const
  {
    const bool needTOF =
      (lstarRequireTOFMatchDe.value != 0) ||
      (lstarEnableBetaCutDe.value != 0) ||
      (lstarEnableExpSignalTOFDe.value != 0 && lstarTOFExpSignalDiffDeMax.value > 0.f);

    if (needTOF) {
      if constexpr (requires { trk.hasTOF(); }) {
        if (!trk.hasTOF()) {
          return false;
        }
      } else {
        return false;
      }
    }

    if (lstarEnableBetaCutDe.value != 0) {
      if constexpr (requires { trk.beta(); }) {
        if (trk.beta() <= lstarBetaCutDe.value) {
          return false;
        }
      } else {
        return false;
      }
    }

    if (lstarEnableExpSignalTOFDe.value != 0 &&
        lstarTOFExpSignalDiffDeMax.value > 0.f) {
      if constexpr (requires { trk.tofExpSignalDiffDe(); }) {
        if (std::abs(trk.tofExpSignalDiffDe()) > lstarTOFExpSignalDiffDeMax.value) {
          return false;
        }
      } else {
        return false;
      }
    }

    return true;
  }

  bool passFinalCandidatePID(float pt, float nsTPC, float nsTOF, bool hasTof,
                             float tpcCut, float tofCut, float circularCut,
                             float ptRef, bool isDeuteron = false) const
  {
    // Strategy 1: analysis-note style circular TPC+TOF cut
    if (lstarPidStrategy.value == PidStrategyCircularTPCAndTOF) {
      if (pt < ptRef) {
        return std::abs(nsTPC) < tpcCut;
      }

      if (!hasTof) {
        return false;
      }

      return std::sqrt(nsTPC * nsTPC + nsTOF * nsTOF) < circularCut;
    }

    if (lstarPidStrategy.value == PidStrategyNucleiDeuteronTPC) {
      // For deuterons, use TPC nσ as the main hard PID selection.
      // If the optional deuteron TOF nσ cut is enabled, apply it only when
      // the track has a valid TOF match. Tracks without TOF are kept with
      // the TPC-only decision, avoiding an additional TOF-matching efficiency loss.
      if (isDeuteron) {
        if (std::abs(nsTPC) >= tpcCut) {
          return false;
        }
        if (lstarEnableTOFNsigmaCutDe.value != 0 && hasTof) {
          return std::abs(nsTOF) < tofCut;
        }
        return true;
      }

      // For kaons/protons, use the Lambda(1520)-like pT-ref circular TPC+TOF logic.
      if (pt < ptRef) {
        return std::abs(nsTPC) < tpcCut;
      }

      if (!hasTof) {
        return false;
      }

      return std::sqrt(nsTPC * nsTPC + nsTOF * nsTOF) < circularCut;
    }

    // Strategy 0: pT-ref dependent rectangular TPC+TOF cut.
    // Below pTref use TPC only; above pTref require TOF and apply both TPC and TOF cuts.
    if (pt < ptRef) {
      return std::abs(nsTPC) < tpcCut;
    }
    if (!hasTof) {
      return false;
    }
    return (std::abs(nsTPC) < tpcCut) && (std::abs(nsTOF) < tofCut);
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
    return RecoDecay::m(std::array{std::array{px1, py1, pz1},
                                   std::array{px2, py2, pz2}},
                        std::array{m1, m2});
  }

  template <bool IsMC, typename TCollision, typename TTracks>
  void processReco(TCollision const& collision, TTracks const& tracks)
  {
    histos.fill(HIST("hEventSelection"), 1.);
    // Use unfiltered collisions here so the first bin includes vertex rejects.
    // Preserve the strict vertex acceptance of collisionZVtxFilter; the MC QA
    // process still uses that framework filter and does not fill this counter.
    if (!(std::abs(collision.posZ()) < lstarCutVertex.value)) {
      return;
    }
    // Event selection (cfgTrigger) -- AO2D only
    if (!keepCollisionAO2D(collision)) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 2.);
    // physics masses (GeV/c^2)
    constexpr double MassProton = o2::constants::physics::MassProton;
    constexpr double MassKaonCharged = o2::constants::physics::MassKaonCharged;

    std::vector<KaonCand> kaonCands;
    std::vector<ProxyCand> proxyCands;
    std::vector<ProtonCand> protonCands;
    kaonCands.reserve(128);
    proxyCands.reserve(32);
    protonCands.reserve(128);

    float eventMultFallback = 0.f; // fallback mixing variable: number of selected charged tracks (after quality cuts)

    // Inclusive track loop before final candidate PID cuts.
    // It counts selected charged tracks for fallback multiplicity and when enabled,
    // fills inclusive PID QA after track-quality cuts but before final K/p/d PID cuts.
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

      const double pForPidQA = static_cast<double>(trk.pt()) * std::cosh(static_cast<double>(trk.eta()));
      const bool hasTofForPidQA = hasTOFMatch(trk);

      histos.fill(HIST("hHasTOFVsP"), pForPidQA, hasTofForPidQA ? 1.0 : 0.0);
      histos.fill(HIST("hHasTOFVsPt"), trk.pt(), hasTofForPidQA ? 1.0 : 0.0);

      // Per-species PID-QA (tagged) histograms
      const int sp = classifyPidSpecies(trk);
      switch (sp) {
        case 0: { // pion
          histos.fill(HIST("hNsigmaTPCPionTaggedVsP"), pForPidQA, trk.tpcNSigmaPi());
          histos.fill(HIST("hNsigmaTOFPionTaggedVsP"), pForPidQA, trk.tofNSigmaPi());
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
          histos.fill(HIST("hNsigmaTPCKaonTaggedVsP"), pForPidQA, trk.tpcNSigmaKa());
          histos.fill(HIST("hNsigmaTOFKaonTaggedVsP"), pForPidQA, trk.tofNSigmaKa());
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
          histos.fill(HIST("hNsigmaTPCProtonTaggedVsP"), pForPidQA, trk.tpcNSigmaPr());
          histos.fill(HIST("hNsigmaTOFProtonTaggedVsP"), pForPidQA, trk.tofNSigmaPr());
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
          histos.fill(HIST("hNsigmaTPCDeuteronTaggedVsP"), pForPidQA, trk.tpcNSigmaDe());
          histos.fill(HIST("hNsigmaTOFDeuteronTaggedVsP"), pForPidQA, trk.tofNSigmaDe());
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
    // Exactly once after event cuts, before PID probes or candidate requirements.
    fillEventMultiplicityCounters(collision, eventMult);
    if (lstarEnableProxyPIDProbes.value != 0) {
      if constexpr (IsMC) {
        fillProxyPIDProbes<true>(tracks, eventMult, collision.mcCollisionId());
      } else {
        fillProxyPIDProbes<false>(tracks, eventMult, -1);
      }
    }

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

      // Deuteron kinematics needed before PID because the PID strategy can depend on pT
      const float ptD = trkD.pt();
      const float etaD = trkD.eta();
      const float phiD = trkD.phi();

      // PID for proxy candidates.
      // Normal mode: use the standard deuteron PID and build p_proxy = p_d / 2.
      // Closure/control mode: select the subset of standard deuteron-proxy candidates
      // that are also proton-like. This tests proton contamination inside the deuteron
      // selection, not all proton candidates.
      const bool useProtonAsProxy = (lstarProxyUseProtonPIDAsDeuteron.value != 0);

      const float nsTPCDe = trkD.tpcNSigmaDe();
      const float nsTOFDe = trkD.tofNSigmaDe();
      // These controls require the actual match flag, never an inferred PID value.
      const bool hasTofDe = trkD.hasTOF();

      if (!passOptionalDeuteronTOFExtras(trkD)) {
        continue;
      }

      const bool passesDeuteronSelection = passFinalCandidatePID(ptD, nsTPCDe, nsTOFDe, hasTofDe,
                                                                 lstarCutNsigmaTPCDe.value, lstarCutNsigmaTOFDe.value,
                                                                 lstarPidCircularCutDe.value, lstarPidPtRefDe.value,
                                                                 true);
      if (!passesDeuteronSelection) {
        continue;
      }

      const float nsTPCPrAsProxy =
        trkD.tpcNSigmaPr();

      const float nsTOFPrAsProxy =
        trkD.tofNSigmaPr();

      const bool passesProtonSelection =
        passFinalCandidatePID(
          ptD,
          nsTPCPrAsProxy,
          nsTOFPrAsProxy,
          hasTofDe,
          lstarCutNsigmaTPCPr.value,
          lstarCutNsigmaTOFPr.value,
          lstarPidCircularCutPr.value,
          lstarPidPtRefPr.value,
          false);

      if (useProtonAsProxy &&
          !passesProtonSelection) {
        continue;
      }

      const double pD = static_cast<double>(ptD) * std::cosh(static_cast<double>(etaD));

      // Candidate QA after final deuteron PID cut
      if (lstarEnablePidQA.value != 0) {
        histos.fill(HIST("hDeuteronProxyPt"), ptD);
        histos.fill(HIST("hDeuteronProxyEta"), etaD);
        histos.fill(HIST("hDeuteronProxyPhi"), phiD);
        histos.fill(HIST("hNsigmaTPCDeuteron"), nsTPCDe);
        histos.fill(HIST("hNsigmaTOFDeuteron"), nsTOFDe);
        histos.fill(HIST("hNsigmaTPCDeuteronVsP"), pD, nsTPCDe);
        histos.fill(HIST("hNsigmaTOFDeuteronVsP"), pD, nsTOFDe);
        histos.fill(HIST("hNsigmaTPCDeuteronVsPt"), ptD, nsTPCDe);
        histos.fill(HIST("hNsigmaTOFDeuteronVsPt"), ptD, nsTOFDe);
        histos.fill(
          HIST(
            "hNsigmaTPCProtonForSelectedDeuteronVsPt"),
          ptD,
          nsTPCPrAsProxy);

        histos.fill(
          HIST(
            "hNsigmaTOFProtonForSelectedDeuteronVsPt"),
          ptD,
          nsTOFPrAsProxy);

        histos.fill(
          HIST(
            "hTPCNsigmaDeVsPrForSelectedDeuteron"),
          nsTPCDe,
          nsTPCPrAsProxy);

        if (hasTofDe) {
          histos.fill(
            HIST(
              "hTOFNsigmaDeVsPrForSelectedDeuteron"),
            nsTOFDe,
            nsTOFPrAsProxy);
          histos.fill(HIST("hTPCvsTOFNsigma_D"), nsTPCDe, nsTOFDe);
        }
        if constexpr (requires { trkD.beta(); }) {
          const float beta = trkD.beta();
          if (hasTofDe && beta > TofBetaMin && beta < TofBetaMax) {
            histos.fill(HIST("hTOFBetaVsPt_D"), ptD, beta);
          }
        }
        if constexpr (requires { trkD.dcaXY(); trkD.dcaZ(); }) {
          histos.fill(HIST("hDCAxyVsPt_D"), ptD, trkD.dcaXY());
          histos.fill(HIST("hDCAzVsPt_D"), ptD, trkD.dcaZ());
        }
      }

      // build proton-proxy momentum from deuteron: p_p ≈ p_d / 2
      const float pxProxy = ProxyMomentumScale * ptD * std::cos(phiD);
      const float pyProxy = ProxyMomentumScale * ptD * std::sin(phiD);
      const float pzProxy = ProxyMomentumScale * ptD * std::sinh(etaD);

      // Full measured momentum of the same d-selected track.
      // This will later be interpreted under the proton mass hypothesis.
      const float pxFull =
        ptD * std::cos(phiD);

      const float pyFull =
        ptD * std::sin(phiD);

      const float pzFull =
        ptD * std::sinh(etaD);

      if (lstarEnableSparse.value != 0 && lstarEnableProxyControls.value != 0) {
        histos.fill(HIST("hLambdaStarProxySelectedTrackControlSparse"), ptD, hasTofDe,
                    passesProtonSelection, std::abs(nsTPCPrAsProxy) < lstarCutNsigmaTPCPr.value);
      }

      proxyCands.push_back(
        ProxyCand{
          .px = pxProxy,
          .py = pyProxy,
          .pz = pzProxy,

          .pxFull = pxFull,
          .pyFull = pyFull,
          .pzFull = pzFull,

          .protonLike =
            passesProtonSelection,
          .hasTOF = hasTofDe,
          .protonTPCCompatible = std::abs(nsTPCPrAsProxy) < lstarCutNsigmaTPCPr.value,

          .charge =
            static_cast<int>(
              trkD.sign()),

          .tid = static_cast<int>(trkD.globalIndex()),
          .truth = {}});
      if constexpr (IsMC) {
        proxyCands.back().truth = trackTruth(trkD, collision.mcCollisionId(), lstarEnableMCProxySpecies.value != 0);
      }
    }

    // Proton candidates (for genuine pK #Lambda^{*} reconstruction)
    for (auto const& trkP : tracks) {
      if (trkP.pt() < lstarCutPtMin.value || std::abs(trkP.eta()) > lstarCutEtaMax.value) {
        continue;
      }
      if (!passTrackQuality(trkP)) {
        continue;
      }
      if (trkP.sign() == 0) {
        continue;
      }

      const float ptP = trkP.pt();
      const float etaP = trkP.eta();
      const float phiP = trkP.phi();

      const float nsTPCPr = trkP.tpcNSigmaPr();
      const float nsTOFPr = trkP.tofNSigmaPr();
      const bool hasTofPr = hasTOFMatch(trkP);
      const bool isProton = passFinalCandidatePID(ptP,
                                                  nsTPCPr,
                                                  nsTOFPr,
                                                  hasTofPr,
                                                  lstarCutNsigmaTPCPr.value,
                                                  lstarCutNsigmaTOFPr.value,
                                                  lstarPidCircularCutPr.value,
                                                  lstarPidPtRefPr.value, false);
      if (!isProton) {
        continue;
      }

      const double pP = static_cast<double>(ptP) * std::cosh(static_cast<double>(etaP));

      if (lstarEnablePidQA.value != 0) {
        histos.fill(HIST("hNsigmaTPCProtonVsP"), pP, nsTPCPr);
        histos.fill(HIST("hNsigmaTOFProtonVsP"), pP, nsTOFPr);
        histos.fill(HIST("hNsigmaTPCProtonVsPt"), ptP, nsTPCPr);
        histos.fill(HIST("hNsigmaTOFProtonVsPt"), ptP, nsTOFPr);
        if (hasTofPr) {
          histos.fill(HIST("hTPCvsTOFNsigma_P"), nsTPCPr, nsTOFPr);
        }
        if constexpr (requires { trkP.beta(); }) {
          const float beta = trkP.beta();
          if (hasTofPr && beta > TofBetaMin && beta < TofBetaMax) {
            histos.fill(HIST("hTOFBetaVsPt_P"), ptP, beta);
          }
        }
        if constexpr (requires { trkP.dcaXY(); trkP.dcaZ(); }) {
          histos.fill(HIST("hDCAxyVsPt_P"), ptP, trkP.dcaXY());
          histos.fill(HIST("hDCAzVsPt_P"), ptP, trkP.dcaZ());
        }
      }

      const float pxP = ptP * std::cos(phiP);
      const float pyP = ptP * std::sin(phiP);
      const float pzP = ptP * std::sinh(etaP);

      protonCands.push_back(ProtonCand{
        .px = pxP,
        .py = pyP,
        .pz = pzP,
        .charge = static_cast<int>(trkP.sign()),
        .tid = static_cast<int>(trkP.globalIndex()),
        .truth = {}});
      if constexpr (IsMC) {
        protonCands.back().truth = trackTruth(trkP, collision.mcCollisionId());
      }
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

      // Kaon kinematics needed before PID because the PID strategy can depend on pT
      const float ptK = trkK.pt();
      const float etaK = trkK.eta();
      const float phiK = trkK.phi();

      // PID for kaon candidates
      const float nsTPCK = trkK.tpcNSigmaKa();
      const float nsTOFK = trkK.tofNSigmaKa();
      const bool hasTofK = hasTOFMatch(trkK);
      const bool isKaon = passFinalCandidatePID(ptK,
                                                nsTPCK,
                                                nsTOFK,
                                                hasTofK,
                                                lstarCutNsigmaTPCKaon.value,
                                                lstarCutNsigmaTOFKaon.value,
                                                lstarPidCircularCutKaon.value,
                                                lstarPidPtRefKaon.value, false);
      if (!isKaon) {
        continue;
      }

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
        histos.fill(HIST("hNsigmaTPCKaonVsPt"), ptK, nsTPCK);
        histos.fill(HIST("hNsigmaTOFKaonVsPt"), ptK, nsTOFK);
        if (hasTofK) {
          histos.fill(HIST("hTPCvsTOFNsigma_K"), nsTPCK, nsTOFK);
        }
        if constexpr (requires { trkK.beta(); }) {
          const float beta = trkK.beta();
          if (hasTofK && beta > TofBetaMin && beta < TofBetaMax) {
            histos.fill(HIST("hTOFBetaVsPt_K"), ptK, beta);
          }
        }
        if constexpr (requires { trkK.dcaXY(); trkK.dcaZ(); }) {
          histos.fill(HIST("hDCAxyVsPt_K"), ptK, trkK.dcaXY());
          histos.fill(HIST("hDCAzVsPt_K"), ptK, trkK.dcaZ());
        }
      }

      const float pxK = ptK * std::cos(phiK);
      const float pyK = ptK * std::sin(phiK);
      const float pzK = ptK * std::sinh(etaK);

      kaonCands.push_back(KaonCand{
        .px = pxK,
        .py = pyK,
        .pz = pzK,
        .hasTOF = hasTofK,
        .charge = static_cast<int>(trkK.sign()),
        .tid = static_cast<int>(trkK.globalIndex()),
        .truth = {}});
      if constexpr (IsMC) {
        kaonCands.back().truth = trackTruth(trkK, collision.mcCollisionId(), lstarEnableMCProxySpecies.value != 0);
      }
    }

    if (kaonCands.empty()) {
      // Still update mixing buffer so that later events can mix with this one.
      // The pK mixed-event background needs stored protons, while the proxy background
      // needs stored proxy candidates. Therefore, do not require proxy candidates here.
      if (mixingEnabled()) {
        LStarMixEventEntry entry;
        entry.mult = eventMult;
        entry.zvtx = collision.posZ();
        entry.kaons = std::move(kaonCands);
        entry.protons = std::move(protonCands);

        mLStarMixEvents.push_front(std::move(entry));
        if (mLStarMixEvents.size() > static_cast<size_t>(lstarNoMixedEvents.value)) {
          mLStarMixEvents.pop_back();
        }
      }
      return;
    }

    const bool hasProtonCandidates = !protonCands.empty();
    const bool hasProxyCandidates = !proxyCands.empty();

    // --- SAME-EVENT: genuine pK #Lambda^{*} candidates ---
    if (hasProtonCandidates) {
      for (auto const& pr : protonCands) {
        for (auto const& k : kaonCands) {
          if (pr.tid == k.tid) {
            continue;
          }
          const double mass = invariantMass(pr.px, pr.py, pr.pz, MassProton,
                                            k.px, k.py, k.pz, MassKaonCharged);

          const float pxTot = pr.px + k.px;
          const float pyTot = pr.py + k.py;
          const float pzTot = pr.pz + k.pz;
          const float ptPair = ptFromPxPy(pxTot, pyTot);

          const float yPair = rapidityFromMomentumAndMass(pxTot, pyTot, pzTot, mass);

          if (std::abs(yPair) > lstarLambdaAbsYMax.value) {
            continue;
          }

          if constexpr (IsMC) {
            fillRecoMC(pr, k, mass, ptPair, false);
          }
          const bool unlikeSignPK = (pr.charge * k.charge) < 0;
          if (unlikeSignPK) {
            histos.fill(HIST("hInvMassPKUnlike"), mass);
            if (lstarEnableSparse.value != 0) {
              histos.fill(HIST("hLambdaStarPKUnlikeSparse"), mass, ptPair, eventMult);
            }
          } else {
            histos.fill(HIST("hInvMassPKLike"), mass);
            if (lstarEnableSparse.value != 0) {
              histos.fill(HIST("hLambdaStarPKLikeSparse"), mass, ptPair, eventMult);
            }
          }
        }
      }
    }

    // --- SAME-EVENT: proxy (d/2) + K ---
    if (hasProxyCandidates) {
      for (auto const& pr : proxyCands) {
        for (auto const& k : kaonCands) {
          if (pr.tid == k.tid) {
            continue; // sanity check: should never match, but just in case of bug in candidate-building logic
          }
          const double mass = invariantMass(pr.px, pr.py, pr.pz, MassProton, k.px, k.py, k.pz, MassKaonCharged);
          // Same d-selected track, but now use its FULL measured
          // momentum and assign the proton mass.
          //
          // No proton PID requirement is involved.
          const double massFullPK =
            invariantMass(
              pr.pxFull,
              pr.pyFull,
              pr.pzFull,
              MassProton,
              k.px,
              k.py,
              k.pz,
              MassKaonCharged);

          const float pxTotFullPK =
            pr.pxFull + k.px;

          const float pyTotFullPK =
            pr.pyFull + k.py;

          const float pzTotFullPK =
            pr.pzFull + k.pz;

          const float yFullPK =
            rapidityFromMomentumAndMass(
              pxTotFullPK,
              pyTotFullPK,
              pzTotFullPK,
              massFullPK);

          const float pxTot = pr.px + k.px;
          const float pyTot = pr.py + k.py;
          const float pzTot = pr.pz + k.pz;
          const float ptPair = ptFromPxPy(pxTot, pyTot);

          const float yProxy =
            rapidityFromMomentumAndMass(
              pxTot,
              pyTot,
              pzTot,
              mass);

          if (std::abs(yProxy) >
              lstarLambdaAbsYMax.value) {
            continue;
          }

          if constexpr (IsMC) {
            fillRecoMC(pr, k, mass, ptPair, true);
            if (lstarEnableMCProxySpecies.value != 0) {
              fillMCProxySpecies(pr, k, mass, ptPair, eventMult);
            }
          }
          const bool unlikeSignProxy =
            (pr.charge * k.charge) < 0;

          // Inclusive invariant-mass spectrum for the #Lambda^{*} proxy (d/2 + K)
          histos.fill(HIST("hDeuteronProxyMass"), mass);
          if (lstarEnableSparse.value != 0) {
            if (lstarEnableProxyControls.value != 0) {
              const bool passFullY = std::abs(yFullPK) <= lstarLambdaAbsYMax.value;
              if (unlikeSignProxy) {
                histos.fill(HIST("hLambdaStarProxyControlUnlikeSparse"), mass, ptPair, eventMult, pr.hasTOF, pr.protonLike, pr.protonTPCCompatible);
                if (lstarEnableFullPKMaps.value != 0) {
                  histos.fill(HIST("hLambdaStarProxyControlFullPKUnlikeSparse"), mass, ptPair, eventMult, pr.hasTOF, pr.protonLike, pr.protonTPCCompatible, massFullPK, passFullY);
                }
              } else {
                histos.fill(HIST("hLambdaStarProxyControlLikeSparse"), mass, ptPair, eventMult, pr.hasTOF, pr.protonLike, pr.protonTPCCompatible);
                if (lstarEnableFullPKMaps.value != 0) {
                  histos.fill(HIST("hLambdaStarProxyControlFullPKLikeSparse"), mass, ptPair, eventMult, pr.hasTOF, pr.protonLike, pr.protonTPCCompatible, massFullPK, passFullY);
                }
              }
            }

            if (unlikeSignProxy) {

              histos.fill(
                HIST(
                  "hLambdaStarProxySparse"),
                mass,
                ptPair,
                eventMult);

            } else {

              histos.fill(
                HIST(
                  "hLambdaStarProxyLikeSparse"),
                mass,
                ptPair,
                eventMult);
            }
            if (std::abs(yFullPK) <=
                lstarLambdaAbsYMax.value) {

              if (unlikeSignProxy) {

                if (lstarEnableFullPKMaps.value != 0) {
                  histos.fill(
                    HIST(
                      "hLambdaStarProxyVsFullPKUnlikeSparse"),
                    mass,
                    massFullPK,
                    ptPair,
                    eventMult);
                }

              } else {

                if (lstarEnableFullPKMaps.value != 0) {
                  histos.fill(
                    HIST(
                      "hLambdaStarProxyVsFullPKLikeSparse"),
                    mass,
                    massFullPK,
                    ptPair,
                    eventMult);
                }
              }
            }
            if (pr.protonLike) {

              if (unlikeSignProxy) {

                if (lstarEnableProtonLikeSparse.value != 0) {
                  histos.fill(
                    HIST("hLambdaStarProxyProtonLikeUnlikeSparse"),
                    mass,
                    ptPair,
                    eventMult);
                }

              } else {

                if (lstarEnableProtonLikeSparse.value != 0) {
                  histos.fill(
                    HIST("hLambdaStarProxyProtonLikeLikeSparse"),
                    mass,
                    ptPair,
                    eventMult);
                }
              }
            }
          }
        }
      }
    }

    if (mixingEnabled()) {
      // --- MIXED-EVENT: current kaons + previous-event real protons ---
      // This fills the standard pK mixed-event background.
      for (auto const& prev : mLStarMixEvents) {

        if (std::abs(prev.zvtx - collision.posZ()) >
            lstarMixZvtxMax.value) {
          continue;
        }

        if (std::abs(prev.mult - eventMult) >
            lstarMixMultMax.value) {
          continue;
        }

        if (prev.protons.empty()) {
          continue;
        }

        for (auto const& pr : prev.protons) {

          for (auto const& k : kaonCands) {

            // Unlike-sign pK mixed-event background.
            if ((pr.charge * k.charge) >= 0) {
              continue;
            }

            const double mass =
              invariantMass(
                pr.px,
                pr.py,
                pr.pz,
                MassProton,
                k.px,
                k.py,
                k.pz,
                MassKaonCharged);

            const float pxTot =
              pr.px + k.px;

            const float pyTot =
              pr.py + k.py;

            const float pzTot =
              pr.pz + k.pz;

            const float ptPair =
              ptFromPxPy(
                pxTot,
                pyTot);

            const float yPair =
              rapidityFromMomentumAndMass(
                pxTot,
                pyTot,
                pzTot,
                mass);

            if (std::abs(yPair) >
                lstarLambdaAbsYMax.value) {
              continue;
            }

            if (lstarEnableSparse.value != 0) {

              if (lstarNoMixedEvents.value > 0) {
                histos.fill(
                  HIST(
                    "hLambdaStarPKMixedSparse"),
                  mass,
                  ptPair,
                  eventMult);
              }
            }
          }
        }
      }

      // --- MIXED-EVENT: current proxies + previous-event kaons ---
      // This fills the deuteron-proxy mixed-event background.
      //
      // For every d-selected proxy candidate we calculate:
      //
      //   1) M_{(d/2)K}
      //      using p_proxy = p_track / 2
      //
      //   2) M_{pK}^{full}
      //      using the FULL measured momentum of the same d-selected track,
      //      but assigning the proton mass.
      //
      // The standard proxy mixed spectrum requires only the proxy rapidity
      // acceptance. The 2D proxy-vs-full-pK map additionally requires the
      // full-pK hypothesis to satisfy the same rapidity acceptance.
      if (hasProxyCandidates) {

        for (auto const& prev : mLStarMixEvents) {

          if (std::abs(prev.zvtx - collision.posZ()) >
              lstarMixZvtxMax.value) {
            continue;
          }

          if (std::abs(prev.mult - eventMult) >
              lstarMixMultMax.value) {
            continue;
          }

          if (prev.kaons.empty()) {
            continue;
          }

          for (auto const& pr : proxyCands) {

            for (auto const& k : prev.kaons) {

              // Unlike-sign proxy-K mixed-event background.
              if ((pr.charge * k.charge) >= 0) {
                continue;
              }

              // ----------------------------------------------------------
              // Proxy hypothesis: p_proxy = p_d / 2
              // ----------------------------------------------------------

              const double mass =
                invariantMass(
                  pr.px,
                  pr.py,
                  pr.pz,
                  MassProton,
                  k.px,
                  k.py,
                  k.pz,
                  MassKaonCharged);

              const float pxTot =
                pr.px + k.px;

              const float pyTot =
                pr.py + k.py;

              const float pzTot =
                pr.pz + k.pz;

              const float ptPair =
                ptFromPxPy(
                  pxTot,
                  pyTot);

              const float yProxy =
                rapidityFromMomentumAndMass(
                  pxTot,
                  pyTot,
                  pzTot,
                  mass);

              // Standard proxy rapidity acceptance.
              if (std::abs(yProxy) >
                  lstarLambdaAbsYMax.value) {
                continue;
              }

              // ----------------------------------------------------------
              // Full-momentum proton hypothesis for the SAME d-selected
              // track.
              // ----------------------------------------------------------

              const double massFullPK =
                invariantMass(
                  pr.pxFull,
                  pr.pyFull,
                  pr.pzFull,
                  MassProton,
                  k.px,
                  k.py,
                  k.pz,
                  MassKaonCharged);

              const float pxTotFullPK =
                pr.pxFull + k.px;

              const float pyTotFullPK =
                pr.pyFull + k.py;

              const float pzTotFullPK =
                pr.pzFull + k.pz;

              const float yFullPK =
                rapidityFromMomentumAndMass(
                  pxTotFullPK,
                  pyTotFullPK,
                  pzTotFullPK,
                  massFullPK);

              if (lstarEnableSparse.value != 0) {

                if (lstarEnableProxyControls.value != 0) {
                  const bool passFullY = std::abs(yFullPK) <= lstarLambdaAbsYMax.value;
                  histos.fill(HIST("hLambdaStarProxyControlMixedSparse"), mass, ptPair, eventMult, pr.hasTOF, pr.protonLike, pr.protonTPCCompatible);
                  if (lstarEnableFullPKMaps.value != 0) {
                    histos.fill(HIST("hLambdaStarProxyControlFullPKMixedSparse"), mass, ptPair, eventMult, pr.hasTOF, pr.protonLike, pr.protonTPCCompatible, massFullPK, passFullY);
                  }
                }

                // Standard mixed-event proxy spectrum.
                //
                // This requires only the proxy hypothesis to satisfy
                // |y| < lstarLambdaAbsYMax.
                if (lstarNoMixedEvents.value > 0) {
                  histos.fill(
                    HIST(
                      "hLambdaStarProxyMixedSparse"),
                    mass,
                    ptPair,
                    eventMult);
                }

                // Mixed-event 2D map.
                //
                // Since this histogram explicitly represents both
                // hypotheses, require the full-pK hypothesis to satisfy
                // the same configurable rapidity acceptance as well.
                if (std::abs(yFullPK) <=
                    lstarLambdaAbsYMax.value) {

                  if (lstarEnableFullPKMaps.value != 0) {
                    if (lstarNoMixedEvents.value > 0) {
                      histos.fill(
                        HIST(
                          "hLambdaStarProxyVsFullPKMixedSparse"),
                        mass,
                        massFullPK,
                        ptPair,
                        eventMult);
                    }
                  }
                }
              }
            }
          }
        }
      }

      // --- Update mixing buffer with current event ---
      LStarMixEventEntry entry;
      entry.mult = eventMult;
      entry.zvtx = collision.posZ();
      entry.kaons = std::move(kaonCands);
      entry.protons = std::move(protonCands);

      mLStarMixEvents.push_front(std::move(entry));
      if (mLStarMixEvents.size() > static_cast<size_t>(lstarNoMixedEvents.value)) {
        mLStarMixEvents.pop_back();
      }
    }
  }

  void processData(CollisionsWithEvSel::iterator const& collision, FilteredTracks const& tracks)
  {
    processReco<false>(collision, tracks);
  }

  // A literal process() is always invoked by O2 in addition to optional methods.
  // Keep the established JSON key "process" while registering data only once.
  PROCESS_SWITCH_FULL(Lambdastarproxy, processData, process, "Lambda* proxy via (d/2)+K", true);

  void processMCReco(CollisionsWithEvSelMC::iterator const& collision,
                     FilteredTracksMC const& tracks, aod::McParticles const&,
                     aod::McCollisions const&)
  {
    // Inclusive reconstructed spectra include unmatched tracks and collisions,
    // just as in data. Truth components require a valid same-event MC label.
    processReco<true>(collision, tracks);
  }
  PROCESS_SWITCH(Lambdastarproxy, processMCReco, "MC: data-identical reconstruction and truth-component spectra", false);

  void processMCGen(aod::McCollision const& collision, aod::McParticles const& particles)
  {
    fillGenerated(collision, particles);
  }
  PROCESS_SWITCH(Lambdastarproxy, processMCGen, "MC: generated pK and (d/2)K spectra, without reconstructed-event requirements", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<Lambdastarproxy>(context)};
}
