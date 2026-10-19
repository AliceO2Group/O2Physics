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
//
/// \file twoParticleCorrelationsMultSpher.cxx
/// \brief Two-particle angular correlations in pp collisions at 13 TeV.
///        Structured in parallel with the AliRoot AliMESppColTask analysis.
///        MC truth level: McCollision + McParticles.
///        Includes: multiplicity classes, transverse sphericity, PID (MC: PDG code),
///                  event mixing with 4 separate pools per sphericity class,
///                  THnSparse for multi-dimensional analysis + TH2D SE/ME for
///                  direct normalisation (both per-trigger and per-pair).
///
/// Pool structure (separate for MC generated and Reco reconstructed data):
///   eventPoolsMC[sphClass][zvtxBin][multBin]
///   eventPoolsReco[sphClass][zvtxBin][multBin]
///   sphClass 0 : S_T in (0.0, 0.3]  — jetty
///   sphClass 1 : S_T in (0.3, 0.6]  — intermediate
///   sphClass 2 : S_T in (0.6, 1.0]  — isotropic
///   sphClass 3 : S_T >  0.0          — all events (no sphericity cut)
///
/// \author Madalina Tarzila

#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
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
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPDGCode.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// ============================================================
//  PID species
// ============================================================
enum PidSpecies { PidUnidentified = 0,
                  PidPion,
                  PidKaon,
                  PidProton,
                  NPidSpecies };
static constexpr std::array<const char*, NPidSpecies> PidNames{"unid", "pion", "kaon", "proton"};

// ============================================================
//  Sphericity classes
// ============================================================
static constexpr int NSphClasses = 4;
static constexpr std::array<double, NSphClasses> SphMax{0.3, 0.6, 1.0, 1.0};
static constexpr std::array<const char*, NSphClasses> SphLabels{"jetty", "intermediate", "isotropic", "all"};

// Half a bin: puts integer-valued axes (indices, cut-flow steps) at bin centres.
static constexpr double HalfBin = 0.5;

// ============================================================
//  Cut-flow steps (bin numbers of the evSel_* and track cut-flow histograms)
// ============================================================
enum EventCut { EventCutRead = 1,
                EventCutZvtx,
                EventCutMultBin,
                EventCutZvtxBin,
                EventCutLeadPt,
                EventCutSphericity,
                EventCutSphClass,
                EventCutSeFilled,
                NEventCuts = EventCutSeFilled };
enum TrackCut { TrackCutSeen = 1,
                TrackCutGlobal,
                TrackCutDcaXY,
                TrackCutDcaZ,
                TrackCutEta,
                TrackCutPt,
                NTrackCuts = TrackCutPt };

// ============================================================
struct TwoParticleCorrelationsMultSpher {

  // ----------------------------------------------------------------
  //  Configurables
  // ----------------------------------------------------------------
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histograms [0-10 GeV/c]"};
  Configurable<int> nBinsPt2{"nBinsPt2", 50, "N bins in pT study histograms [0-5 GeV/c]"};
  Configurable<int> mixPoolSize{"mixPoolSize", 200, "Max events per mixing pool cell"};
  Configurable<double> minAssocPt{"minAssocPt", 0.2, "Min pT for track selection (GeV/c)"};
  Configurable<double> minLeadPt{"minLeadPt", 1.0, "Min pT for leading (trigger) particle (GeV/c)"};
  Configurable<double> maxLeadPt{"maxLeadPt", 2.0, "Max pT for leading (trigger) particle (GeV/c)"};
  Configurable<double> minAssocCorr{"minAssocCorr", 1.0, "Min pT for associated particle in correlations (GeV/c)"};
  Configurable<double> maxAssocCorr{"maxAssocCorr", 2.0, "Max pT for associated particle in correlations (GeV/c)"};
  Configurable<double> maxEta{"maxEta", 0.8, "Max |eta| for track selection"};
  Configurable<double> maxZvtx{"maxZvtx", 10.0, "Max |z_vtx| (cm)"};
  // |y| cut applied ONLY to PID-identified tracks (pion/kaon/proton) --
  // independent of maxEta, which applies to all-species tracks (hSE/hME).
  // Rapidity is well-defined only for massive/identified particles; at
  // pT=0.5 GeV/c, proton y_max~0.48, so |y|<0.5 keeps full acceptance for
  // that species. Reference: ALICE 2511.10399 (arXiv, Nov 2025), Sec. 4 --
  // "tracks reconstructed within |eta|<0.8, with the additional rapidity
  // cut of |y|<0.5".
  Configurable<double> maxRapPID{"maxRapPID", 0.5, "Max |y| for PID-identified tracks (pion/kaon/proton), independent of maxEta"};

  // Reco-specific configurables
  Configurable<double> maxDcaXY{"maxDcaXY", 2.0, "Max DCA XY (cm) for reco tracks"};
  Configurable<double> maxDcaZ{"maxDcaZ", 2.0, "Max DCA Z (cm) for reco tracks"};
  // Combined TPC+TOF PID criterion per ALICE 2511.10399 Sec. 4: two
  // different thresholds -- selection at 2 sigma, ambiguity rejection
  // at 3 sigma.
  Configurable<double> maxNSigmaPID{"maxNSigmaPID", 2.0, "Max combined N_sigma,PID for species selection (ALICE 2511.10399)"};
  Configurable<double> maxNSigmaPIDAmbig{"maxNSigmaPIDAmbig", 3.0, "N_sigma,PID threshold for ambiguity rejection (>1 species below this -> PidUnidentified)"};

  // Thresholds of the legacy PID gating (independent 3-sigma cuts on TPC
  // and TOF, no combined quadrature, no pT-dependent TOF gating), kept so
  // that it can be switched on with useLegacyNSigmaGating without removing
  // the ALICE 2511.10399 criterion.
  Configurable<double> maxNSigmaTPC{"maxNSigmaTPC", 3.0, "[legacy] Max |N_sigma| TPC for PID, used only when useLegacyNSigmaGating=true"};
  Configurable<double> maxNSigmaTOF{"maxNSigmaTOF", 3.0, "[legacy] Max |N_sigma| TOF for PID, used only when useLegacyNSigmaGating=true"};

  // Three independently toggleable flags to compare the legacy PID
  // methodology (independent 3-sigma TPC/TOF cuts, minimum-chi2 species
  // choice, no rapidity cut) with the ALICE 2511.10399-based one. Each
  // isolates one component, so any combination of the 8 can be run to
  // see which change drives a statistics difference. The gating and
  // selection defaults reproduce the 2511.10399 behaviour; setting all
  // three to the legacy side reproduces the legacy pidFromNSigma.
  //
  // applyRapidityCutPID defaults to false: the |y|<0.5 cut only discarded
  // already-correctly-identified tracks between |y|=0.5 and the |eta|<0.8
  // acceptance, so dropping it raises both yield and purity at the same
  // time. Gating and selection stay at the 2511.10399 values, because
  // relaxing the gating trades purity for statistics.
  Configurable<bool> useLegacyNSigmaGating{"useLegacyNSigmaGating", false, "false=pT-gated combined N_sigma,PID (2511.10399); true=independent TPC/TOF 3-sigma AND-gate, TOF always required when available (legacy)"};
  Configurable<bool> useLegacyPIDSelection{"useLegacyPIDSelection", false, "false=reject if >1 species below ambiguity threshold; true=pick minimum-chi2 species among those passing the gate, no rejection (legacy)"};
  Service<o2::framework::O2DatabasePDG> pdg{};

  Configurable<bool> applyRapidityCutPID{"applyRapidityCutPID", false, "false=skip the |y|<maxRapPID cut entirely for PID-identified tracks (the legacy methodology had no such cut); true=apply it"};

  // ----------------------------------------------------------------
  //  Per-track container stored in mixing pool
  //  NOTE: no PID in pool — ME does not use PID
  // ----------------------------------------------------------------
  struct TrackSimple {
    float phi = 0.f;
    float eta = 0.f;
    float pt = 0.f;
    int pdg = 0; // PDG code: used only for SE sparse/TH2; 0 for reco
  };

  // ----------------------------------------------------------------
  //  Event class constants
  // ----------------------------------------------------------------
  static constexpr float MaxValidNSigma = 100.f;  // larger |nSigma| means the detector had no signal (-999)
  static constexpr float MinPtCombinedPid = 0.5f; // GeV/c, above this TPC+TOF are combined
  static constexpr std::size_t MinTracksSphericity = 3;
  static constexpr double Epsilon = 1e-12;
  static constexpr double MinCharge = 0.1; // charged particle: |q| in units of e
  static constexpr int NMultBins = 6;
  static constexpr int NZvtxBins = 5;
  // Bin edges: bin i spans [edge[i], edge[i+1]); values outside [front, back) give -1 in findBin.
  static constexpr std::array<double, NZvtxBins + 1> ZvtxEdges{-10., -5., -2.5, 2.5, 5., 10.};
  static constexpr std::array<int, NMultBins + 1> MultEdges{1, 10, 15, 20, 30, 46, 81};

  static std::string multLabel(int im)
  {
    return std::to_string(MultEdges[im]) + "-" + std::to_string(MultEdges[im + 1] - 1);
  }

  // ----------------------------------------------------------------
  //  Mixing pools
  //  [sphClass][zvtxBin][multBin]     — "all" (existing)
  //  [pidIdx][sphClass][zvtxBin][multBin] — per PID species
  //    pidIdx: 0=pion, 1=kaon, 2=proton
  // ----------------------------------------------------------------
  static constexpr int NPidCorr = 3; // pion, kaon, proton
  static constexpr std::array<const char*, NPidCorr> PidCorrNames{"pion", "kaon", "proton"};
  // PDG codes corresponding to pidIdx
  static int pidCorrPDG(int idx)
  {
    if (idx == 0) {
      return PDG_t::kPiPlus;
    }
    if (idx == 1) {
      return PDG_t::kKPlus;
    }
    return PDG_t::kProton;
  }
  /// Index of a pion/kaon/proton in the NPidCorr-sized arrays, or -1 for any other species.
  static int pidSpeciesToIdx(PidSpecies pid)
  {
    if (pid < PidPion || pid > PidProton) {
      return -1;
    }
    return pid - PidPion;
  }

  using EventPool = std::deque<std::vector<TrackSimple>>;
  // Pool "all" — separate per MC (generated) and Reco (reconstructed)
  std::array<std::array<std::array<EventPool, NMultBins>, NZvtxBins>, NSphClasses> eventPoolsMC;
  std::array<std::array<std::array<EventPool, NMultBins>, NZvtxBins>, NSphClasses> eventPoolsReco;
  // Pool per PID species (like-sign) — separate for MC and Reco
  std::array<std::array<std::array<std::array<EventPool, NMultBins>, NZvtxBins>, NSphClasses>, NPidCorr> eventPoolsPidMC;
  std::array<std::array<std::array<std::array<EventPool, NMultBins>, NZvtxBins>, NSphClasses>, NPidCorr> eventPoolsPidReco;

  // ----------------------------------------------------------------
  //  Histogram registry + handles
  // ----------------------------------------------------------------
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // SE/ME TH2D "all": [sphClass][multBin] — MC and Reco separate
  std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses> hSeMC{};
  std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses> hMeMC{};
  std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses> hSeReco{};
  std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses> hMeReco{};

  // SE/ME TH2D like-sign per species: [pidIdx][sphClass][multBin] — MC and Reco separate
  std::array<std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses>, NPidCorr> hSePidMC{};
  std::array<std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses>, NPidCorr> hMePidMC{};
  std::array<std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses>, NPidCorr> hSePidReco{};
  std::array<std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses>, NPidCorr> hMePidReco{};

  // Trigger counter: axes = mult bin (x) x sphClass (y) — MC and Reco separate
  std::shared_ptr<TH2> hNtrigMC;
  std::shared_ptr<TH2> hNtrigReco;
  std::array<std::shared_ptr<TH2>, NPidCorr> hNtrigPidMC{};
  std::array<std::shared_ptr<TH2>, NPidCorr> hNtrigPidReco{};

  // Per-mult sphericity — MC and Reco separate
  std::array<std::shared_ptr<TH1>, NMultBins> hSphMultMC{};
  std::array<std::shared_ptr<TH1>, NMultBins> hSphMultReco{};

  // ----------------------------------------------------------------
  //  [DEBUG, TEMPORARY] Track-cut counter for processReco diagnostics.
  //  Bins: 1=tracks seen, 2=passed isGlobalTrackWoDCA, 3=passed DCAxy,
  //        4=passed DCAz, 5=passed eta, 6=passed pt (selected).
  //  Remove once the empty-Reco-histograms issue is understood/fixed.
  // ----------------------------------------------------------------
  std::shared_ptr<TH1> hTrackCutDebugReco;

  // ----------------------------------------------------------------
  //  PID confusion matrix (true species from PDG vs. tagged species from
  //  pidFromNSigma applied to the SAME reconstructed track, per pT bin).
  //  Quick contamination estimate, reported as a systematic; it does NOT
  //  correct anything. Filled only by processMCConfusion (needs both MC
  //  truth and reconstructed nSigma together, unlike
  //  processMC/processReco which only ever see one or the other).
  // ----------------------------------------------------------------
  std::shared_ptr<TH3> hPidConfusionMC;

  // ----------------------------------------------------------------
  //  Helpers
  // ----------------------------------------------------------------

  static PidSpecies pidFromPDG(int pdg)
  {
    switch (std::abs(pdg)) {
      case PDG_t::kPiPlus:
        return PidPion;
      case PDG_t::kKPlus:
        return PidKaon;
      case PDG_t::kProton:
        return PidProton;
      default:
        return PidUnidentified;
    }
  }

  /// PID from TPC(+TOF) nSigma for reconstructed tracks. Two independently
  /// toggleable behaviors (see the useLegacyNSigmaGating and
  /// useLegacyPIDSelection Configurables), so that the legacy and the
  /// ALICE 2511.10399-based methodologies can be compared:
  ///
  /// legacyGating (gate = "does this species hypothesis pass at all"):
  ///   false (default): pT>0.5 && TOF available -> combined
  ///     N_sigma,PID = sqrt(nsTPC^2+nsTOF^2) < maxSel (2 sigma, per ALICE
  ///     2511.10399 Sec. 4); else TPC-only.
  ///   true: legacy style -- |nsTPC| < maxTPC AND (if TOF available) |nsTOF| < maxTOF,
  ///     independently (no quadrature sum for the gate itself), TOF
  ///     required whenever available regardless of pT.
  ///
  /// legacySelection (what to do with hypotheses that pass the gate):
  ///   false (default): reject the track outright (PidUnidentified) if MORE
  ///     THAN ONE hypothesis also falls below the looser ambiguity
  ///     threshold maxAmbig (3 sigma) -- this comparison always uses the
  ///     combined-quadrature N_sigma,PID, regardless of legacyGating, so
  ///     the two flags stay independently meaningful.
  ///   true: legacy style -- among hypotheses passing the gate,
  ///     pick the one with minimum chi2 = nsTPC^2+nsTOF^2 (or nsTPC^2
  ///     alone if no TOF); no rejection for multiple passing hypotheses.
  ///
  /// The two flags combine independently (e.g. legacy gating + new
  /// selection is a valid, meaningful hybrid) -- this is deliberate, so
  /// the statistics impact of each component can be isolated rather than
  /// only comparing "old vs new" as one indivisible block.
  template <typename T>
  PidSpecies pidFromNSigma(const T& track, double maxSel, double maxAmbig,
                           double maxTPC, double maxTOF,
                           bool legacyGating, bool legacySelection) const
  {
    struct Hypo {
      PidSpecies pid = PidUnidentified;
      double nsTPC = 0.;
      double nsTOF = 0.;
    };
    const std::array<Hypo, NPidCorr> hypos{{{PidPion, track.tpcNSigmaPi(), track.tofNSigmaPi()},
                                            {PidKaon, track.tpcNSigmaKa(), track.tofNSigmaKa()},
                                            {PidProton, track.tpcNSigmaPr(), track.tofNSigmaPr()}}};

    // TOF unavailable -> nSigma returns -999; check with any species
    bool hasTOF = (std::fabs(track.tofNSigmaPi()) < MaxValidNSigma);
    bool useCombined = (track.pt() > MinPtCombinedPid) && hasTOF; // only feeds the new-style gate

    auto passesGate = [&](const Hypo& h) {
      if (legacyGating) {
        return (std::fabs(h.nsTPC) < maxTPC) && (!hasTOF || std::fabs(h.nsTOF) < maxTOF);
      }
      double nsPID = useCombined
                       ? std::sqrt(h.nsTPC * h.nsTPC + h.nsTOF * h.nsTOF)
                       : std::fabs(h.nsTPC);
      return nsPID < maxSel;
    };

    if (legacySelection) {
      // Best-fit (min chi2) among hypotheses passing the gate -- no
      // ambiguity rejection (legacy style).
      PidSpecies best = PidUnidentified;
      double bestChi2 = std::numeric_limits<double>::max();
      for (const auto& h : hypos) {
        if (!passesGate(h)) {
          continue;
        }
        double chi2 = h.nsTPC * h.nsTPC + (hasTOF ? h.nsTOF * h.nsTOF : 0.0);
        if (chi2 < bestChi2) {
          bestChi2 = chi2;
          best = h.pid;
        }
      }
      return best;
    }

    // Ambiguity-rejection selection (current default style).
    int nBelowAmbig = 0;
    PidSpecies selected = PidUnidentified;
    for (const auto& h : hypos) {
      if (!passesGate(h)) {
        continue;
      }
      double nsPID = useCombined
                       ? std::sqrt(h.nsTPC * h.nsTPC + h.nsTOF * h.nsTOF)
                       : std::fabs(h.nsTPC);
      if (nsPID < maxAmbig) {
        ++nBelowAmbig;
      }
      selected = h.pid;
    }
    if (nBelowAmbig > 1) {
      return PidUnidentified; // ambiguous, reject
    }
    return selected;
  }

  /// Returns PDG mass for identified species, 0 for unidentified.
  static double pidMass(PidSpecies pid)
  {
    switch (pid) {
      case PidPion:
        return o2::constants::physics::MassPiPlus;
      case PidKaon:
        return o2::constants::physics::MassKPlus;
      case PidProton:
        return o2::constants::physics::MassProton;
      default:
        return 0.;
    }
  }

  // Returns specific sph class (0,1,2) or -1 if unphysical
  // Note: class 3 ("all") is handled separately — it is always filled
  // when sph > 0, without needing a specific bin lookup.
  static int sphClass(double sph)
  {
    if (sph <= 0.) {
      return -1;
    }
    for (int ic = 0; ic < NSphClasses - 1; ++ic) {
      if (sph <= SphMax[ic]) {
        return ic;
      }
    }
    return -1;
  }

  static double computeSphericity(const std::vector<TrackSimple>& tracks)
  {
    if (tracks.size() < MinTracksSphericity) {
      return -1.0;
    }
    double sxx = 0., sxy = 0., syy = 0.;
    for (const auto& t : tracks) {
      double px = t.pt * std::cos(t.phi);
      double py = t.pt * std::sin(t.phi);
      sxx += px * px;
      sxy += px * py;
      syy += py * py;
    }
    double tr = sxx + syy;
    double disc = tr * tr - 4.0 * (sxx * syy - sxy * sxy);
    if (disc < 0.) {
      disc = 0.;
    }
    double lambdaMin = (tr - std::sqrt(disc)) / 2.0;
    if (tr < Epsilon) {
      return -1.0;
    }
    return std::clamp(2.0 * lambdaMin / tr, 0.0, 1.0);
  }

  /// Compute rapidity from pT, eta and particle mass.
  /// y = 0.5 * ln((energy + pz) / (energy - pz))
  /// For unidentified particles (mass <= 0), returns eta.
  static double computeRapidity(double pt, double eta, double mass)
  {
    if (mass <= 0.) {
      return eta; // massless / unidentified -> y = eta
    }
    double pz = pt * std::sinh(eta);
    double energy = std::sqrt(mass * mass + pt * pt + pz * pz);
    if (energy <= std::fabs(pz) + Epsilon) {
      return eta; // safety
    }
    return 0.5 * std::log((energy + pz) / (energy - pz));
  }

  // ----------------------------------------------------------------
  //  init
  // ----------------------------------------------------------------
  void init(InitContext const&)
  {
    const AxisSpec axisZvtx{100, -15., 15., "z_{vtx} (cm)"};
    const AxisSpec axisMult{100, 0., 180., "N_{ch}"};
    const AxisSpec axisEta{30, -1.5, 1.5, "#eta"};
    const AxisSpec axisPhi{60, 0., o2::constants::math::TwoPI, "#phi"};
    const AxisSpec axisDphi{30, -o2::constants::math::PIHalf, o2::constants::math::PI + o2::constants::math::PIHalf, "#Delta#phi"};
    const AxisSpec axisDeta{32, -1.6, 1.6, "#Delta#eta"};
    // Reserved for a possible rapidity-differentiated analysis -- unused
    // while the PID correlations run pseudorapidity-only.
    [[maybe_unused]] const AxisSpec axisDy{32, -1.6, 1.6, "#Delta y"};
    const AxisSpec axisPt{nBinsPt, 0., 10., "p_{T} (GeV/c)"};
    const AxisSpec axisSph{100, 0., 1., "S_{T}"};
    const AxisSpec axisPt2{nBinsPt2, 0., 5., "p_{T} (GeV/c)"};

    // ----------------------------------------------------------------
    //  Common lambda: books a complete set of histograms for one
    //  "channel" (MC or Reco), with all names suffixed accordingly.
    //  Takes references to all channel-specific members.
    // ----------------------------------------------------------------
    auto bookChannel = [&](const char* suf,
                           std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses>& hSE,
                           std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses>& hME,
                           std::array<std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses>, NPidCorr>& hSEpid,
                           std::array<std::array<std::array<std::shared_ptr<TH2>, NMultBins>, NSphClasses>, NPidCorr>& hMEpid,
                           std::shared_ptr<TH2>& hNtrig,
                           std::array<std::shared_ptr<TH2>, NPidCorr>& hNtrigPID,
                           std::array<std::shared_ptr<TH1>, NMultBins>& hSphMult) {
      // QA
      auto hQA = histos.add<TH1>(Form("evSel_%s", suf), Form("Event selection (%s)", suf), HistType::kTH1D, {{NEventCuts, HalfBin, static_cast<double>(NEventCuts) + HalfBin}});
      hQA->GetXaxis()->SetBinLabel(EventCutRead, "Events read");
      hQA->GetXaxis()->SetBinLabel(EventCutZvtx, "|z_{vtx}|<10");
      hQA->GetXaxis()->SetBinLabel(EventCutMultBin, "multBin OK");
      hQA->GetXaxis()->SetBinLabel(EventCutZvtxBin, "zvtxBin OK");
      hQA->GetXaxis()->SetBinLabel(EventCutLeadPt, "leadPt #in [min,max]");
      hQA->GetXaxis()->SetBinLabel(EventCutSphericity, "S_{T} valid");
      hQA->GetXaxis()->SetBinLabel(EventCutSphClass, "sphClass OK");
      hQA->GetXaxis()->SetBinLabel(EventCutSeFilled, "SE filled");

      histos.add(Form("Zvertex_%s", suf), Form("z_{vtx} (%s);z_{vtx} (cm);events", suf), kTH1D, {axisZvtx});
      histos.add(Form("multHist_%s", suf), Form("N_{ch} (%s);N_{ch};events", suf), kTH1D, {axisMult});
      histos.add(Form("etaHistogram_%s", suf), Form("#eta (%s);#eta;tracks", suf), kTH1D, {axisEta});
      histos.add(Form("phiHistogram_%s", suf), Form("#phi (%s);#phi;tracks", suf), kTH1D, {axisPhi});
      histos.add(Form("ptHistogram_%s", suf), Form("p_{T} (%s);p_{T} (GeV/c);tracks", suf), kTH1D, {axisPt});
      histos.add(Form("ptLeadHistogram_%s", suf), Form("p_{T} lead (%s);p_{T} (GeV/c);events", suf), kTH1D, {axisPt});
      histos.add(Form("ptAssocHistogram_%s", suf), Form("p_{T} assoc (%s);p_{T} (GeV/c);pairs", suf), kTH1D, {axisPt});
      histos.add(Form("sphericity_%s", suf), Form("S_{T} (inclusive, %s);S_{T};events", suf), kTH1D, {axisSph});
      histos.add(Form("sphericity_beforeCuts_%s", suf), Form("S_{T} (all evts with >=3 tracks, %s);S_{T};events", suf), kTH1D, {axisSph});

      // Per-mult sphericity
      for (int im = 0; im < NMultBins; ++im) {
        hSphMult[im] = histos.add<TH1>(
          Form("sphericity_mult%d_%s", im, suf),
          Form("S_T [%s] (%s);S_{T};events", multLabel(im).c_str(), suf),
          HistType::kTH1D, {axisSph});
      }

      // Trigger counter TH2: x = multBin, y = sphClass
      hNtrig = histos.add<TH2>(
        Form("hNtrig%s", suf),
        Form("N_{trig} (%s);mult bin;sph class", suf),
        HistType::kTH2D,
        {{NMultBins, -HalfBin, static_cast<double>(NMultBins) - HalfBin},
         {NSphClasses, -HalfBin, static_cast<double>(NSphClasses) - HalfBin}});
      for (int im = 0; im < NMultBins; ++im) {
        hNtrig->GetXaxis()->SetBinLabel(im + 1, multLabel(im).c_str());
      }
      for (int ic = 0; ic < NSphClasses; ++ic) {
        hNtrig->GetYaxis()->SetBinLabel(ic + 1, SphLabels[ic]);
      }

      // SE/ME TH2D: [sphClass][multBin]
      for (int ic = 0; ic < NSphClasses; ++ic) {
        for (int im = 0; im < NMultBins; ++im) {
          hSE[ic][im] = histos.add<TH2>(
            Form("hSE%s_sph%d_mult%d", suf, ic, im),
            Form("SE [%s][%s] (%s);#Delta#phi;#Delta#eta", SphLabels[ic], multLabel(im).c_str(), suf),
            HistType::kTH2D, {axisDphi, axisDeta});
          hME[ic][im] = histos.add<TH2>(
            Form("hME%s_sph%d_mult%d", suf, ic, im),
            Form("ME [%s][%s] (%s);#Delta#phi;#Delta#eta", SphLabels[ic], multLabel(im).c_str(), suf),
            HistType::kTH2D, {axisDphi, axisDeta});
        }
      }

      // ---- SE/ME like-sign per PID species: [pidIdx][sphClass][multBin] ---
      // PID correlations use Delta-eta, same as the inclusive hSE/hME
      // histograms, NOT Delta-y. The rapidity infrastructure
      // (computeRapidity(), the per-track y field, axisDy) is kept for a
      // possible rapidity-differentiated analysis and is deliberately
      // unused for now.
      for (int ip = 0; ip < NPidCorr; ++ip) {
        for (int ic = 0; ic < NSphClasses; ++ic) {
          for (int im = 0; im < NMultBins; ++im) {
            hSEpid[ip][ic][im] = histos.add<TH2>(
              Form("hSE%s_pid%d_sph%d_mult%d", suf, ip, ic, im),
              Form("SE %s--%s [%s][%s] (%s);#Delta#varphi;#Delta#eta",
                   PidCorrNames[ip], PidCorrNames[ip], SphLabels[ic], multLabel(im).c_str(), suf),
              HistType::kTH2D, {axisDphi, axisDeta});
            hMEpid[ip][ic][im] = histos.add<TH2>(
              Form("hME%s_pid%d_sph%d_mult%d", suf, ip, ic, im),
              Form("ME %s--%s [%s][%s] (%s);#Delta#varphi;#Delta#eta",
                   PidCorrNames[ip], PidCorrNames[ip], SphLabels[ic], multLabel(im).c_str(), suf),
              HistType::kTH2D, {axisDphi, axisDeta});
          }
        }
        // Trigger counter per PID species
        hNtrigPID[ip] = histos.add<TH2>(
          Form("hNtrig%s_pid%d", suf, ip),
          Form("N_{trig} %s (%s);mult bin;sph class", PidCorrNames[ip], suf),
          HistType::kTH2D,
          {{NMultBins, -HalfBin, static_cast<double>(NMultBins) - HalfBin},
           {NSphClasses, -HalfBin, static_cast<double>(NSphClasses) - HalfBin}});
        for (int im = 0; im < NMultBins; ++im) {
          hNtrigPID[ip]->GetXaxis()->SetBinLabel(im + 1, multLabel(im).c_str());
        }
        for (int ic = 0; ic < NSphClasses; ++ic) {
          hNtrigPID[ip]->GetYaxis()->SetBinLabel(ic + 1, SphLabels[ic]);
        }
      }
    };

    // Book histograms for MC (generated) channel
    bookChannel("MC",
                hSeMC, hMeMC, hSePidMC, hMePidMC,
                hNtrigMC, hNtrigPidMC, hSphMultMC);

    // Book histograms for Reco (reconstructed) channel
    bookChannel("Reco",
                hSeReco, hMeReco, hSePidReco, hMePidReco,
                hNtrigReco, hNtrigPidReco, hSphMultReco);

    // [DEBUG, TEMPORARY] track-cut diagnostic counter for processReco
    hTrackCutDebugReco = histos.add<TH1>(
      "hTrackCutDebug_Reco",
      "Track cut flow (Reco, DEBUG);cut;tracks",
      HistType::kTH1D, {{NTrackCuts, HalfBin, static_cast<double>(NTrackCuts) + HalfBin}});
    hTrackCutDebugReco->GetXaxis()->SetBinLabel(TrackCutSeen, "seen");
    hTrackCutDebugReco->GetXaxis()->SetBinLabel(TrackCutGlobal, "passed isGlobalTrackWoDCA");
    hTrackCutDebugReco->GetXaxis()->SetBinLabel(TrackCutDcaXY, "passed DCAxy");
    hTrackCutDebugReco->GetXaxis()->SetBinLabel(TrackCutDcaZ, "passed DCAz");
    hTrackCutDebugReco->GetXaxis()->SetBinLabel(TrackCutEta, "passed eta");
    hTrackCutDebugReco->GetXaxis()->SetBinLabel(TrackCutPt, "passed pt (selected)");

    // PID confusion matrix (contamination estimate) --
    // axes: pT, true species (PDG), tagged species (pidFromNSigma).
    // Species axis order matches PidSpecies: 0=unid,1=pion,2=kaon,3=proton.
    const AxisSpec axisPidSpecies{static_cast<int>(NPidSpecies), -HalfBin, static_cast<double>(NPidSpecies) - HalfBin, "species"};
    hPidConfusionMC = histos.add<TH3>(
      "hPidConfusion_MC",
      "True vs. tagged PID species (MC);p_{T} (GeV/c);true species;tagged species",
      HistType::kTH3D, {axisPt2, axisPidSpecies, axisPidSpecies});
    for (int ip = 0; ip < static_cast<int>(NPidSpecies); ++ip) {
      hPidConfusionMC->GetYaxis()->SetBinLabel(ip + 1, PidNames[ip]);
      hPidConfusionMC->GetZaxis()->SetBinLabel(ip + 1, PidNames[ip]);
    }
  }

  // ----------------------------------------------------------------
  //  processMC
  // ----------------------------------------------------------------
  void processMC(aod::McCollision const& mcCollision,
                 aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("evSel_MC"), EventCutRead);

    const double zvtx = mcCollision.posZ();
    if (std::fabs(zvtx) > maxZvtx) {
      return;
    }
    histos.fill(HIST("evSel_MC"), EventCutZvtx);
    histos.fill(HIST("Zvertex_MC"), zvtx);

    // ---- one-pass: build selTracks + find leading ---------------
    std::vector<TrackSimple> selTracks;
    selTracks.reserve(MultEdges.back() - 1);
    int leadIdx = -1;
    double pTlead = -1., phiLead = 0., etaLead = 0.;

    for (const auto& p : mcParticles) {
      if (!p.isPhysicalPrimary()) {
        continue;
      }
      auto* pdgPtr = pdg->GetParticle(p.pdgCode());
      if (!pdgPtr || std::fabs(pdgPtr->Charge()) < MinCharge) {
        continue;
      }
      if (std::fabs(p.eta()) > maxEta) {
        continue;
      }
      if (p.pt() < minAssocPt) {
        continue;
      }

      double mass = pdgPtr->Mass();
      double rap = computeRapidity(p.pt(), p.eta(), mass);

      // |y| < maxRapPID cut applies ONLY to PID-identified species
      // (pion/kaon/proton) -- the all-species channel (hSE/hME) keeps the
      // track regardless, gated only by |eta| < maxEta above. Clean way to
      // do this: keep pdg for all-species bookkeeping, but zero it (->
      // PidUnidentified via pidFromPDG) when the rapidity cut fails, so the
      // track contributes to hSE/hME but not hSEpid/hMEpid.
      int pdgForPID = p.pdgCode();
      if (applyRapidityCutPID && pidFromPDG(pdgForPID) != PidUnidentified && std::fabs(rap) > maxRapPID) {
        pdgForPID = 0;
      }

      int idx = static_cast<int>(selTracks.size());
      selTracks.push_back({static_cast<float>(p.phi()), static_cast<float>(p.eta()), static_cast<float>(p.pt()), pdgForPID});

      histos.fill(HIST("etaHistogram_MC"), p.eta());
      histos.fill(HIST("phiHistogram_MC"), p.phi());
      histos.fill(HIST("ptHistogram_MC"), p.pt());

      if (p.pt() > pTlead) {
        pTlead = p.pt();
        phiLead = p.phi();
        etaLead = p.eta();
        leadIdx = idx;
      }
    }

    const int nch = static_cast<int>(selTracks.size());
    histos.fill(HIST("multHist_MC"), nch);

    const int mBin = o2::analysis::findBin(&MultEdges, nch);
    if (mBin < 0) {
      return;
    }
    histos.fill(HIST("evSel_MC"), EventCutMultBin);

    const int zBin = o2::analysis::findBin(&ZvtxEdges, zvtx);
    if (zBin < 0) {
      return;
    }
    histos.fill(HIST("evSel_MC"), EventCutZvtxBin);

    // ---- leading particle pT cut [minLeadPt, maxLeadPt] -----------
    if (leadIdx < 0) {
      return;
    }
    if (pTlead < minLeadPt || pTlead > maxLeadPt) {
      return;
    }
    histos.fill(HIST("evSel_MC"), EventCutLeadPt);
    histos.fill(HIST("ptLeadHistogram_MC"), pTlead);

    // ---- sphericity ---------------------------------------------
    const double sph = computeSphericity(selTracks);
    // Fill sphericity_beforeCuts here — no additional cuts applied yet
    if (sph >= 0.) {
      histos.fill(HIST("sphericity_beforeCuts_MC"), sph);
    }
    if (sph < 0.) {
      return; // fewer than 3 tracks
    }
    histos.fill(HIST("evSel_MC"), EventCutSphericity);
    histos.fill(HIST("sphericity_MC"), sph);
    hSphMultMC[mBin]->Fill(sph);

    const int sphCls = sphClass(sph); // 0,1,2 or -1
    if (sphCls < 0) {
      return;
    }
    histos.fill(HIST("evSel_MC"), EventCutSphClass);

    const PidSpecies pidLead = pidFromPDG(selTracks[leadIdx].pdg);

    // ---- trigger counter ----------------------------------------
    hNtrigMC->Fill(static_cast<double>(mBin), static_cast<double>(sphCls));
    hNtrigMC->Fill(static_cast<double>(mBin), static_cast<double>((NSphClasses - 1))); // "all"

    // ================================================================
    //  SAME-EVENT correlations
    // ================================================================
    bool seHasPairs = false;
    for (int ia = 0; ia < nch; ++ia) {
      if (ia == leadIdx) {
        continue;
      }
      const auto& assoc = selTracks[ia];
      if (assoc.pt >= pTlead) {
        continue; // associated pT < leading pT
      }
      if (assoc.pt < minAssocCorr || assoc.pt > maxAssocCorr) {
        continue;
      }

      histos.fill(HIST("ptAssocHistogram_MC"), assoc.pt);

      const PidSpecies pidAssoc = pidFromPDG(assoc.pdg);

      const double dphi = RecoDecay::constrainAngle(phiLead - assoc.phi, -o2::constants::math::PIHalf);
      const double deta = etaLead - assoc.eta;

      // TH2D SE: specific sph class (uses Delta-eta)
      hSeMC[sphCls][mBin]->Fill(dphi, deta);
      // TH2D SE: "all" class (uses Delta-eta)
      hSeMC[NSphClasses - 1][mBin]->Fill(dphi, deta);

      // SE like-sign: fill only if trigger and associated are the same species
      // Uses Delta-eta (pseudorapidity-only regime)
      int pidLeadIdx = pidSpeciesToIdx(pidLead);
      int pidAssocIdx = pidSpeciesToIdx(pidAssoc);
      if (pidLeadIdx >= 0 && pidLeadIdx == pidAssocIdx) {
        hSePidMC[pidLeadIdx][sphCls][mBin]->Fill(dphi, deta);
        hSePidMC[pidLeadIdx][NSphClasses - 1][mBin]->Fill(dphi, deta);
      }

      seHasPairs = true;
    }
    if (seHasPairs) {
      histos.fill(HIST("evSel_MC"), EventCutSeFilled);
    }

    // ---- hNtrigMC per PID species leading --------------------------
    int pidLeadIdx = pidSpeciesToIdx(pidLead);
    if (pidLeadIdx >= 0) {
      hNtrigPidMC[pidLeadIdx]->Fill(static_cast<double>(mBin), static_cast<double>(sphCls));
      hNtrigPidMC[pidLeadIdx]->Fill(static_cast<double>(mBin), static_cast<double>((NSphClasses - 1)));
    }

    // ================================================================
    //  Update pools AFTER SE
    // ================================================================
    // Pool "all" — specific sph
    {
      auto& pool = eventPoolsMC[sphCls][zBin][mBin];
      pool.push_back(selTracks);
      if (static_cast<int>(pool.size()) > mixPoolSize) {
        pool.pop_front();
      }
    }
    // Pool "all" — all sph
    {
      auto& poolAll = eventPoolsMC[NSphClasses - 1][zBin][mBin];
      poolAll.push_back(selTracks);
      if (static_cast<int>(poolAll.size()) > mixPoolSize) {
        poolAll.pop_front();
      }
    }

    // Pool per PID species — contains only the tracks of that species
    for (int ip = 0; ip < NPidCorr; ++ip) {
      // Build a vector with only the tracks of species ip
      std::vector<TrackSimple> pidTracks;
      for (const auto& t : selTracks) {
        if (pidSpeciesToIdx(pidFromPDG(t.pdg)) == ip) {
          pidTracks.push_back(t);
        }
      }
      if (pidTracks.empty()) {
        continue;
      }
      // Specific sph pool
      {
        auto& pool = eventPoolsPidMC[ip][sphCls][zBin][mBin];
        pool.push_back(pidTracks);
        if (static_cast<int>(pool.size()) > mixPoolSize) {
          pool.pop_front();
        }
      }
      // "all" sph pool
      {
        auto& pool = eventPoolsPidMC[ip][NSphClasses - 1][zBin][mBin];
        pool.push_back(pidTracks);
        if (static_cast<int>(pool.size()) > mixPoolSize) {
          pool.pop_front();
        }
      }
    }

    // ================================================================
    //  MIXED-EVENT "all"
    // ================================================================
    auto doMixing = [&](int poolIdx) {
      const auto& pool = eventPoolsMC[poolIdx][zBin][mBin];
      const int nEvts = static_cast<int>(pool.size()) - 1;
      for (int ie = 0; ie < nEvts; ++ie) {
        for (const auto& tr : pool[ie]) {
          if (tr.pt >= pTlead) {
            continue;
          }
          if (tr.pt < minAssocCorr || tr.pt > maxAssocCorr) {
            continue;
          }
          const double dphi = RecoDecay::constrainAngle(phiLead - tr.phi, -o2::constants::math::PIHalf);
          const double deta = etaLead - tr.eta;
          hMeMC[poolIdx][mBin]->Fill(dphi, deta);
        }
      }
    };

    doMixing(sphCls);
    doMixing(NSphClasses - 1);

    // ================================================================
    //  MIXED-EVENT like-sign per PID species
    //  Trigger = leading track of species pidLeadIdx
    //  Associated = tracks of the same species from the PID pool
    //  Uses Delta-eta (pseudorapidity-only regime)
    // ================================================================
    if (pidLeadIdx >= 0) {
      auto doMixingPID = [&](int poolIdx) {
        const auto& pool = eventPoolsPidMC[pidLeadIdx][poolIdx][zBin][mBin];
        const int nEvts = static_cast<int>(pool.size()) - 1;
        for (int ie = 0; ie < nEvts; ++ie) {
          for (const auto& tr : pool[ie]) {
            if (tr.pt >= pTlead) {
              continue;
            }
            if (tr.pt < minAssocCorr || tr.pt > maxAssocCorr) {
              continue;
            }
            const double dphi = RecoDecay::constrainAngle(phiLead - tr.phi, -o2::constants::math::PIHalf);
            const double detaPID = etaLead - tr.eta;
            hMePidMC[pidLeadIdx][poolIdx][mBin]->Fill(dphi, detaPID);
          }
        }
      };
      doMixingPID(sphCls);
      doMixingPID(NSphClasses - 1);
    }
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMultSpher, processMC, "Process MC truth", true);

  // ================================================================
  //  processReco — reconstructed data (Collision + Tracks)
  // ================================================================
  using ColReco = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using TrackReco = soa::Join<aod::TracksIU, aod::TracksDCA, aod::TrackSelection,
                              aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                              aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;

  void processReco(ColReco::iterator const& col,
                   TrackReco const& tracks)
  {
    histos.fill(HIST("evSel_Reco"), EventCutRead);

    // ---- event selection ----
    if (!col.sel8()) {
      return;
    }
    const double zvtx = col.posZ();
    if (std::fabs(zvtx) > maxZvtx) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), EventCutZvtx);
    histos.fill(HIST("Zvertex_Reco"), zvtx);

    // ---- one-pass: build selTracks + find leading ---------------
    std::vector<TrackSimple> selTracks;
    selTracks.reserve(MultEdges.back() - 1);
    int leadIdx = -1;
    double pTlead = -1., phiLead = 0., etaLead = 0.;

    for (const auto& track : tracks) {
      hTrackCutDebugReco->Fill(TrackCutSeen); // seen

      // Quality selection: global track without DCA (apply DCA separately)
      if (!track.isGlobalTrackWoDCA()) {
        continue;
      }
      hTrackCutDebugReco->Fill(TrackCutGlobal); // passed isGlobalTrackWoDCA

      // DCA cuts
      if (std::fabs(track.dcaXY()) > maxDcaXY) {
        continue;
      }
      hTrackCutDebugReco->Fill(TrackCutDcaXY); // passed DCAxy
      if (std::fabs(track.dcaZ()) > maxDcaZ) {
        continue;
      }
      hTrackCutDebugReco->Fill(TrackCutDcaZ); // passed DCAz

      // Kinematic cuts (same as MC)
      if (std::fabs(track.eta()) > maxEta) {
        continue;
      }
      hTrackCutDebugReco->Fill(TrackCutEta); // passed eta
      if (track.pt() < minAssocPt) {
        continue;
      }
      hTrackCutDebugReco->Fill(TrackCutPt); // passed pt (selected)

      // PID from nSigma (combined TPC+TOF, ALICE 2511.10399, unless the
      // legacy flags below switch to legacy gating/selection)
      PidSpecies pid = pidFromNSigma(track, maxNSigmaPID, maxNSigmaPIDAmbig,
                                     maxNSigmaTPC, maxNSigmaTOF,
                                     useLegacyNSigmaGating, useLegacyPIDSelection);
      double mass = pidMass(pid);
      double rap = computeRapidity(track.pt(), track.eta(), mass);

      // |y| < maxRapPID cut applies ONLY to PID-identified species -- see
      // matching comment in processMC. Track stays in selTracks either way
      // (all-species channel unaffected); only the PID assignment is reset.
      if (applyRapidityCutPID && pid != PidUnidentified && std::fabs(rap) > maxRapPID) {
        pid = PidUnidentified;
      }

      // pdg = 0 for reco (no MC truth); store sign for future like/unlike-sign
      int idx = static_cast<int>(selTracks.size());
      selTracks.push_back({static_cast<float>(track.phi()), static_cast<float>(track.eta()), static_cast<float>(track.pt()), 0});
      // Store PID info in pdg field as signed pidIdx for later use:
      //   sign(track.sign()) * (211 for pion, 321 for kaon, 2212 for proton, 0 for unid)
      const int pidIdx = pidSpeciesToIdx(pid);
      const int pdgLike = pidIdx >= 0 ? pidCorrPDG(pidIdx) : 0;
      selTracks.back().pdg = track.sign() * pdgLike;

      histos.fill(HIST("etaHistogram_Reco"), track.eta());
      histos.fill(HIST("phiHistogram_Reco"), track.phi());
      histos.fill(HIST("ptHistogram_Reco"), track.pt());

      if (track.pt() > pTlead) {
        pTlead = track.pt();
        phiLead = track.phi();
        etaLead = track.eta();
        leadIdx = idx;
      }
    }

    const int nch = static_cast<int>(selTracks.size());
    histos.fill(HIST("multHist_Reco"), nch);

    const int mBin = o2::analysis::findBin(&MultEdges, nch);
    if (mBin < 0) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), EventCutMultBin);

    const int zBin = o2::analysis::findBin(&ZvtxEdges, zvtx);
    if (zBin < 0) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), EventCutZvtxBin);

    // ---- leading particle pT cut [minLeadPt, maxLeadPt] -----------
    if (leadIdx < 0) {
      return;
    }
    if (pTlead < minLeadPt || pTlead > maxLeadPt) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), EventCutLeadPt);
    histos.fill(HIST("ptLeadHistogram_Reco"), pTlead);

    // ---- sphericity ---------------------------------------------
    const double sph = computeSphericity(selTracks);
    if (sph >= 0.) {
      histos.fill(HIST("sphericity_beforeCuts_Reco"), sph);
    }
    if (sph < 0.) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), EventCutSphericity);
    histos.fill(HIST("sphericity_Reco"), sph);
    hSphMultReco[mBin]->Fill(sph);

    const int sphCls = sphClass(sph);
    if (sphCls < 0) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), EventCutSphClass);

    const PidSpecies pidLead = pidFromPDG(selTracks[leadIdx].pdg);

    // ---- trigger counter ----------------------------------------
    hNtrigReco->Fill(static_cast<double>(mBin), static_cast<double>(sphCls));
    hNtrigReco->Fill(static_cast<double>(mBin), static_cast<double>((NSphClasses - 1)));

    // ================================================================
    //  SAME-EVENT correlations (identical logic to processMC)
    // ================================================================
    bool seHasPairs = false;
    for (int ia = 0; ia < nch; ++ia) {
      if (ia == leadIdx) {
        continue;
      }
      const auto& assoc = selTracks[ia];
      if (assoc.pt >= pTlead) {
        continue;
      }
      if (assoc.pt < minAssocCorr || assoc.pt > maxAssocCorr) {
        continue;
      }

      histos.fill(HIST("ptAssocHistogram_Reco"), assoc.pt);

      const PidSpecies pidAssoc = pidFromPDG(assoc.pdg);

      const double dphi = RecoDecay::constrainAngle(phiLead - assoc.phi, -o2::constants::math::PIHalf);
      const double deta = etaLead - assoc.eta;

      // All-species SE (Delta-eta)
      hSeReco[sphCls][mBin]->Fill(dphi, deta);
      hSeReco[NSphClasses - 1][mBin]->Fill(dphi, deta);

      // PID SE — same species, Delta-eta (pseudorapidity-only regime)
      int pidLeadIdx = pidSpeciesToIdx(pidLead);
      int pidAssocIdx = pidSpeciesToIdx(pidAssoc);
      if (pidLeadIdx >= 0 && pidLeadIdx == pidAssocIdx) {
        hSePidReco[pidLeadIdx][sphCls][mBin]->Fill(dphi, deta);
        hSePidReco[pidLeadIdx][NSphClasses - 1][mBin]->Fill(dphi, deta);
      }

      seHasPairs = true;
    }
    if (seHasPairs) {
      histos.fill(HIST("evSel_Reco"), EventCutSeFilled);
    }

    // ---- hNtrigReco per PID species leading --------------------------
    int pidLeadIdx = pidSpeciesToIdx(pidLead);
    if (pidLeadIdx >= 0) {
      hNtrigPidReco[pidLeadIdx]->Fill(static_cast<double>(mBin), static_cast<double>(sphCls));
      hNtrigPidReco[pidLeadIdx]->Fill(static_cast<double>(mBin), static_cast<double>((NSphClasses - 1)));
    }

    // ================================================================
    //  Update pools AFTER SE (identical to processMC)
    // ================================================================
    {
      auto& pool = eventPoolsReco[sphCls][zBin][mBin];
      pool.push_back(selTracks);
      if (static_cast<int>(pool.size()) > mixPoolSize) {
        pool.pop_front();
      }
    }
    {
      auto& poolAll = eventPoolsReco[NSphClasses - 1][zBin][mBin];
      poolAll.push_back(selTracks);
      if (static_cast<int>(poolAll.size()) > mixPoolSize) {
        poolAll.pop_front();
      }
    }

    for (int ip = 0; ip < NPidCorr; ++ip) {
      std::vector<TrackSimple> pidTracks;
      for (const auto& t : selTracks) {
        if (pidSpeciesToIdx(pidFromPDG(t.pdg)) == ip) {
          pidTracks.push_back(t);
        }
      }
      if (pidTracks.empty()) {
        continue;
      }
      {
        auto& pool = eventPoolsPidReco[ip][sphCls][zBin][mBin];
        pool.push_back(pidTracks);
        if (static_cast<int>(pool.size()) > mixPoolSize) {
          pool.pop_front();
        }
      }
      {
        auto& pool = eventPoolsPidReco[ip][NSphClasses - 1][zBin][mBin];
        pool.push_back(pidTracks);
        if (static_cast<int>(pool.size()) > mixPoolSize) {
          pool.pop_front();
        }
      }
    }

    // ================================================================
    //  MIXED-EVENT "all" (Delta-eta)
    // ================================================================
    auto doMixing = [&](int poolIdx) {
      const auto& pool = eventPoolsReco[poolIdx][zBin][mBin];
      const int nEvts = static_cast<int>(pool.size()) - 1;
      for (int ie = 0; ie < nEvts; ++ie) {
        for (const auto& tr : pool[ie]) {
          if (tr.pt >= pTlead) {
            continue;
          }
          if (tr.pt < minAssocCorr || tr.pt > maxAssocCorr) {
            continue;
          }
          const double dphi = RecoDecay::constrainAngle(phiLead - tr.phi, -o2::constants::math::PIHalf);
          const double deta = etaLead - tr.eta;
          hMeReco[poolIdx][mBin]->Fill(dphi, deta);
        }
      }
    };

    doMixing(sphCls);
    doMixing(NSphClasses - 1);

    // ================================================================
    //  MIXED-EVENT PID (Delta-eta, pseudorapidity-only regime)
    // ================================================================
    if (pidLeadIdx >= 0) {
      auto doMixingPID = [&](int poolIdx) {
        const auto& pool = eventPoolsPidReco[pidLeadIdx][poolIdx][zBin][mBin];
        const int nEvts = static_cast<int>(pool.size()) - 1;
        for (int ie = 0; ie < nEvts; ++ie) {
          for (const auto& tr : pool[ie]) {
            if (tr.pt >= pTlead) {
              continue;
            }
            if (tr.pt < minAssocCorr || tr.pt > maxAssocCorr) {
              continue;
            }
            const double dphi = RecoDecay::constrainAngle(phiLead - tr.phi, -o2::constants::math::PIHalf);
            const double detaPID = etaLead - tr.eta;
            hMePidReco[pidLeadIdx][poolIdx][mBin]->Fill(dphi, detaPID);
          }
        }
      };
      doMixingPID(sphCls);
      doMixingPID(NSphClasses - 1);
    }
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMultSpher, processReco, "Process reconstructed data", false);

  // ----------------------------------------------------------------
  //  processMCConfusion -- PID contamination estimate. Fills hPidConfusionMC (pT, true species, tagged
  //  species) on MC only -- needs BOTH the truth PDG (via McTrackLabels)
  //  AND the reconstructed track's nSigma values on the SAME track,
  //  which neither processMC (truth only) nor processReco (no truth
  //  link) has access to, hence a dedicated process function. Reports
  //  contamination as a systematic; does NOT correct anything (no
  //  unfolding is implemented). Off by default, like the other MC-only diagnostics.
  // ----------------------------------------------------------------
  using TrackRecoMC = soa::Join<TrackReco, aod::McTrackLabels>;

  void processMCConfusion(ColReco::iterator const& col,
                          TrackRecoMC const& tracks,
                          aod::McParticles const&)
  {
    if (!col.sel8()) {
      return;
    }
    if (std::fabs(col.posZ()) > maxZvtx) {
      return;
    }

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      auto mcPart = track.mcParticle();
      if (!mcPart.isPhysicalPrimary()) {
        continue;
      }
      auto* pdgPtr = pdg->GetParticle(mcPart.pdgCode());
      if (!pdgPtr || std::fabs(pdgPtr->Charge()) < MinCharge) {
        continue;
      }

      // Same track-level cuts as processReco, so the matrix reflects the
      // actual analyzed track population, not an unfiltered sample.
      if (std::fabs(track.dcaXY()) > maxDcaXY) {
        continue;
      }
      if (std::fabs(track.dcaZ()) > maxDcaZ) {
        continue;
      }
      if (std::fabs(track.eta()) > maxEta) {
        continue;
      }
      if (track.pt() < minAssocPt) {
        continue;
      }

      PidSpecies truePid = pidFromPDG(mcPart.pdgCode());
      PidSpecies taggedPid = pidFromNSigma(track, maxNSigmaPID, maxNSigmaPIDAmbig,
                                           maxNSigmaTPC, maxNSigmaTOF,
                                           useLegacyNSigmaGating, useLegacyPIDSelection);

      // Same rapidity cut as processReco applies to its "tagged" PID --
      // without this, applyRapidityCutPID=false would silently have zero
      // effect on this matrix even though it changes the processReco yields.
      if (applyRapidityCutPID && taggedPid != PidUnidentified) {
        double tagMass = pidMass(taggedPid);
        double tagRap = computeRapidity(track.pt(), track.eta(), tagMass);
        if (std::fabs(tagRap) > maxRapPID) {
          taggedPid = PidUnidentified;
        }
      }

      histos.fill(HIST("hPidConfusion_MC"), track.pt(), static_cast<double>(truePid), static_cast<double>(taggedPid));
    }
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMultSpher, processMCConfusion, "Rapid PID contamination estimate: true (PDG) vs tagged (nSigma) species confusion matrix, MC only", false);

}; // end struct

// ----------------------------------------------------------------
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TwoParticleCorrelationsMultSpher>(cfgc)};
}
