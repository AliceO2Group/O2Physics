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
/// \brief Two-particle angular correlations in pp collisions at 13 TeV.
///        Structured in parallel with the AliRoot AliMESppColTask analysis.
///        MC truth level: McCollision + McParticles.
///        Includes: multiplicity classes, transverse sphericity, PID (MC: PDG code),
///                  event mixing with 4 separate pools per sphericity class,
///                  THnSparse for multi-dimensional analysis + TH2D SE/ME for
///                  direct normalisation (both per-trigger and per-pair).
///
/// Pool structure (separate for MC generated and Reco reconstructed data):
///   eventPools_MC[sphClass][zvtxBin][multBin]
///   eventPools_Reco[sphClass][zvtxBin][multBin]
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
#include <TMath.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// ============================================================
//  PID species
// ============================================================
enum PidSpecies { kPidUnidentified = 0,
                  kPidPion,
                  kPidKaon,
                  kPidProton,
                  kNPidSpecies };
static constexpr std::array<const char*, kNPidSpecies> kPidNames{"unid", "pion", "kaon", "proton"};

// ============================================================
//  Sphericity classes
// ============================================================
static constexpr int nSphClasses = 4;
static constexpr std::array<double, nSphClasses> kSphMax{0.3, 0.6, 1.0, 1.0};
static constexpr std::array<const char*, nSphClasses> kSphLabels{"jetty", "intermediate", "isotropic", "all"};

// ============================================================
struct myExampleTask {

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
  // cut of |y|<0.5" (verified 2026-08-11, quoted directly from the paper).
  Configurable<double> maxRapPID{"maxRapPID", 0.5, "Max |y| for PID-identified tracks (pion/kaon/proton), independent of maxEta"};

  // Reco-specific configurables
  Configurable<double> maxDcaXY{"maxDcaXY", 2.0, "Max DCA XY (cm) for reco tracks"};
  Configurable<double> maxDcaZ{"maxDcaZ", 2.0, "Max DCA Z (cm) for reco tracks"};
  // Combined TPC+TOF PID criterion per ALICE 2511.10399 Sec. 4: two
  // DIFFERENT thresholds, not one -- selection at 2 sigma, ambiguity
  // rejection at 3 sigma (verified directly against the paper text
  // 2026-08-11; an earlier draft of this instruction conflated the two
  // into a single 2-sigma threshold, corrected here).
  Configurable<double> maxNSigmaPID{"maxNSigmaPID", 2.0, "Max combined N_sigma,PID for species selection (ALICE 2511.10399)"};
  Configurable<double> maxNSigmaPIDAmbig{"maxNSigmaPIDAmbig", 3.0, "N_sigma,PID threshold for ambiguity rejection (>1 species below this -> kPidUnidentified)"};

  // Pre-2026-08-11 PID thresholds, restored (not reused for anything else)
  // so the legacy gating style below can be toggled back on without
  // deleting the ALICE 2511.10399 criterion. Values match the original
  // (myExampleTask_BeforeOnlyProtons.cxx): independent 3-sigma cuts on
  // TPC and TOF, no combined quadrature, no pT-dependent TOF gating.
  Configurable<double> maxNSigmaTPC{"maxNSigmaTPC", 3.0, "[legacy] Max |N_sigma| TPC for PID, used only when useLegacyNSigmaGating=true"};
  Configurable<double> maxNSigmaTOF{"maxNSigmaTOF", 3.0, "[legacy] Max |N_sigma| TOF for PID, used only when useLegacyNSigmaGating=true"};

  // Three independently toggleable flags to compare pre-2026-08-11 PID
  // methodology against the current ALICE 2511.10399-based one, for the
  // PID-methodology-vs-statistics study (requested 2026-08-13). Each
  // isolates one component instead of a single "old vs new" switch, so
  // any combination of the 8 can be run to see which change drives a
  // statistics difference. Gating/selection defaults still reproduce the
  // 2511.10399 behavior exactly; setting all three to the legacy side
  // reproduces myExampleTask_BeforeOnlyProtons.cxx's pidFromNSigma exactly.
  //
  // applyRapidityCutPID defaults to false as of 2026-09-03, per the
  // completed PID-methodology-vs-statistics study (49-file MC, see
  // ~/EventShape-Wiki/Physics/PID-Methodology-Statistics-Impact-Study.md):
  // the |y|<0.5 cut only discarded already-correctly-identified tracks
  // between |y|=0.5 and the |eta|<0.8 acceptance, so dropping it raises
  // both yield and purity simultaneously -- no tradeoff. Gating/selection
  // are kept at the current (non-legacy) values by the same decision:
  // gating relaxation was found to be a real purity-for-statistics
  // tradeoff, not a free win, and was declined.
  Configurable<bool> useLegacyNSigmaGating{"useLegacyNSigmaGating", false, "false=pT-gated combined N_sigma,PID (2511.10399); true=independent TPC/TOF 3-sigma AND-gate, TOF always required when available (pre-2026-08-11)"};
  Configurable<bool> useLegacyPIDSelection{"useLegacyPIDSelection", false, "false=reject if >1 species below ambiguity threshold; true=pick minimum-chi2 species among those passing the gate, no rejection (pre-2026-08-11)"};
  Service<o2::framework::O2DatabasePDG> pdg{};

  Configurable<bool> applyRapidityCutPID{"applyRapidityCutPID", false, "false=skip the |y|<maxRapPID cut entirely for PID-identified tracks (pre-2026-08-11 had no such cut) -- default as of 2026-09-03, see PID-methodology study"};

  // ----------------------------------------------------------------
  //  Per-track container stored in mixing pool
  //  NOTE: no PID in pool — ME does not use PID
  // ----------------------------------------------------------------
  struct TrackSimple {
    float phi;
    float eta;
    float pt;
    int pdg; // PDG code: used only for SE sparse/TH2; 0 for reco
  };

  // ----------------------------------------------------------------
  //  Event class constants
  // ----------------------------------------------------------------
  static constexpr int nMultBins = 6;
  static constexpr int nZvtxBins = 5;
  // Bin edges: bin i spans [edge[i], edge[i+1]); values outside [front, back) give -1 in findBin.
  static constexpr std::array<double, nZvtxBins + 1> kZvtxEdges{-10., -5., -2.5, 2.5, 5., 10.};
  static constexpr std::array<int, nMultBins + 1> kMultEdges{1, 10, 15, 20, 30, 46, 81};

  static std::string multLabel(int im)
  {
    return std::to_string(kMultEdges[im]) + "-" + std::to_string(kMultEdges[im + 1] - 1);
  }

  // ----------------------------------------------------------------
  //  Mixing pools
  //  [sphClass][zvtxBin][multBin]     — "all" (esistente)
  //  [pidIdx][sphClass][zvtxBin][multBin] — per specie PID
  //    pidIdx: 0=pion, 1=kaon, 2=proton
  // ----------------------------------------------------------------
  static constexpr int nPidCorr = 3; // pion, kaon, proton
  static constexpr std::array<const char*, nPidCorr> kPidCorrNames{"pion", "kaon", "proton"};
  // PDG codes corrispondenti a pidIdx
  static int pidCorrPDG(int idx)
  {
    if (idx == 0) {
      return 211;
    }
    if (idx == 1) {
      return 321;
    }
    return 2212;
  }
  // Mappa PidSpecies -> pidIdx (-1 se non e' pion/kaon/proton)
  static int pidSpeciesToIdx(PidSpecies pid)
  {
    if (pid == kPidPion) {
      return 0;
    }
    if (pid == kPidKaon) {
      return 1;
    }
    if (pid == kPidProton) {
      return 2;
    }
    return -1;
  }

  using EventPool = std::deque<std::vector<TrackSimple>>;
  // Pool "all" — separate per MC (generated) and Reco (reconstructed)
  std::array<std::array<std::array<EventPool, nMultBins>, nZvtxBins>, nSphClasses> eventPools_MC;
  std::array<std::array<std::array<EventPool, nMultBins>, nZvtxBins>, nSphClasses> eventPools_Reco;
  // Pool per specie PID (like-sign) — separate per MC e Reco
  std::array<std::array<std::array<std::array<EventPool, nMultBins>, nZvtxBins>, nSphClasses>, nPidCorr> eventPoolsPID_MC;
  std::array<std::array<std::array<std::array<EventPool, nMultBins>, nZvtxBins>, nSphClasses>, nPidCorr> eventPoolsPID_Reco;

  // ----------------------------------------------------------------
  //  Histogram registry + handles
  // ----------------------------------------------------------------
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // SE/ME TH2D "all": [sphClass][multBin] — MC e Reco separati
  std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses> hSE_MC{};
  std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses> hME_MC{};
  std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses> hSE_Reco{};
  std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses> hME_Reco{};

  // SE/ME TH2D like-sign per specie: [pidIdx][sphClass][multBin] — MC e Reco separati
  std::array<std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses>, nPidCorr> hSEpid_MC{};
  std::array<std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses>, nPidCorr> hMEpid_MC{};
  std::array<std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses>, nPidCorr> hSEpid_Reco{};
  std::array<std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses>, nPidCorr> hMEpid_Reco{};

  // Trigger counter: axes = mult bin (x) x sphClass (y) — MC e Reco separati
  std::shared_ptr<TH2> hNtrig_MC;
  std::shared_ptr<TH2> hNtrig_Reco;
  std::array<std::shared_ptr<TH2>, nPidCorr> hNtrigPID_MC{};
  std::array<std::shared_ptr<TH2>, nPidCorr> hNtrigPID_Reco{};

  // Per-mult sphericity — MC e Reco separati
  std::array<std::shared_ptr<TH1>, nMultBins> hSphMult_MC{};
  std::array<std::shared_ptr<TH1>, nMultBins> hSphMult_Reco{};

  // ----------------------------------------------------------------
  //  [DEBUG, TEMPORARY] Track-cut counter for processReco diagnostics.
  //  Bins: 1=tracks seen, 2=passed isGlobalTrackWoDCA, 3=passed DCAxy,
  //        4=passed DCAz, 5=passed eta, 6=passed pt (selected).
  //  Remove once the empty-Reco-histograms issue is understood/fixed.
  // ----------------------------------------------------------------
  std::shared_ptr<TH1> hTrackCutDebug_Reco;

  // ----------------------------------------------------------------
  //  PID confusion matrix (true species from PDG vs. tagged species from
  //  pidFromNSigma applied to the SAME reconstructed track, per pT bin).
  //  "Rapid" contamination estimate requested 2026-08-14 -- reports
  //  contamination as a systematic, does NOT correct anything. Filled
  //  only by processMCConfusion (needs both MC truth and reconstructed
  //  nSigma together, unlike processMC/processReco which only ever see
  //  one or the other). See auto-memory / STATUS.md for the "rapid vs
  //  riguros" design discussion this came out of.
  // ----------------------------------------------------------------
  std::shared_ptr<TH3> hPidConfusion_MC;

  // ----------------------------------------------------------------
  //  Helpers
  // ----------------------------------------------------------------

  static PidSpecies pidFromPDG(int pdg)
  {
    switch (std::abs(pdg)) {
      case 211:
        return kPidPion;
      case 321:
        return kPidKaon;
      case 2212:
        return kPidProton;
      default:
        return kPidUnidentified;
    }
  }

  /// PID from TPC(+TOF) nSigma for reconstructed tracks. Two independently
  /// toggleable behaviors (2026-08-13, added to compare pre-2026-08-11 PID
  /// methodology against the current ALICE 2511.10399-based one, without
  /// deleting either -- see useLegacyNSigmaGating/useLegacyPIDSelection
  /// Configurables):
  ///
  /// legacyGating (gate = "does this species hypothesis pass at all"):
  ///   false (default): pT>0.5 && TOF available -> combined
  ///     N_sigma,PID = sqrt(nsTPC^2+nsTOF^2) < maxSel (2 sigma, per ALICE
  ///     2511.10399 Sec. 4); else TPC-only.
  ///   true: pre-2026-08-11 style, matches myExampleTask_BeforeOnlyProtons.cxx
  ///     exactly -- |nsTPC| < maxTPC AND (if TOF available) |nsTOF| < maxTOF,
  ///     independently (no quadrature sum for the gate itself), TOF
  ///     required whenever available regardless of pT.
  ///
  /// legacySelection (what to do with hypotheses that pass the gate):
  ///   false (default): reject the track outright (kPidUnidentified) if MORE
  ///     THAN ONE hypothesis also falls below the looser ambiguity
  ///     threshold maxAmbig (3 sigma) -- this comparison always uses the
  ///     combined-quadrature N_sigma,PID, regardless of legacyGating, so
  ///     the two flags stay independently meaningful.
  ///   true: pre-2026-08-11 style -- among hypotheses passing the gate,
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
      PidSpecies pid;
      double nsTPC;
      double nsTOF;
    };
    const std::array<Hypo, 3> hypos{{{kPidPion, track.tpcNSigmaPi(), track.tofNSigmaPi()},
                                     {kPidKaon, track.tpcNSigmaKa(), track.tofNSigmaKa()},
                                     {kPidProton, track.tpcNSigmaPr(), track.tofNSigmaPr()}}};

    // TOF unavailable -> nSigma returns -999; check with any species
    bool hasTOF = (std::fabs(track.tofNSigmaPi()) < 100.f);
    bool useCombined = (track.pt() > 0.5) && hasTOF; // only feeds the new-style gate

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
      // ambiguity rejection (pre-2026-08-11 style).
      PidSpecies best = kPidUnidentified;
      double bestChi2 = 1e9;
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
    PidSpecies selected = kPidUnidentified;
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
      return kPidUnidentified; // ambiguous, reject
    }
    return selected;
  }

  /// Returns PDG mass for identified species, 0 for unidentified.
  static double pidMass(PidSpecies pid)
  {
    switch (pid) {
      case kPidPion:
        return 0.13957;
      case kPidKaon:
        return 0.49368;
      case kPidProton:
        return 0.93827;
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
    for (int ic = 0; ic < nSphClasses - 1; ++ic) {
      if (sph <= kSphMax[ic]) {
        return ic;
      }
    }
    return -1;
  }

  static double computeSphericity(const std::vector<TrackSimple>& tracks)
  {
    if (tracks.size() < 3) {
      return -1.0;
    }
    double Sxx = 0., Sxy = 0., Syy = 0.;
    for (const auto& t : tracks) {
      double px = t.pt * std::cos(t.phi);
      double py = t.pt * std::sin(t.phi);
      Sxx += px * px;
      Sxy += px * py;
      Syy += py * py;
    }
    double tr = Sxx + Syy;
    double disc = tr * tr - 4.0 * (Sxx * Syy - Sxy * Sxy);
    if (disc < 0.) {
      disc = 0.;
    }
    double lambdaMin = (tr - std::sqrt(disc)) / 2.0;
    if (tr < 1e-12) {
      return -1.0;
    }
    return std::clamp(2.0 * lambdaMin / tr, 0.0, 1.0);
  }

  /// Compute rapidity from pT, eta and particle mass.
  /// y = 0.5 * ln((E + pz) / (E - pz))
  /// For unidentified particles (mass <= 0), returns eta.
  static double computeRapidity(double pt, double eta, double mass)
  {
    if (mass <= 0.) {
      return eta; // massless / unidentified -> y = eta
    }
    double pz = pt * std::sinh(eta);
    double E = std::sqrt(mass * mass + pt * pt + pz * pz);
    if (E <= std::fabs(pz) + 1e-12) {
      return eta; // safety
    }
    return 0.5 * std::log((E + pz) / (E - pz));
  }

  // ----------------------------------------------------------------
  //  init
  // ----------------------------------------------------------------
  void init(InitContext const&)
  {
    const AxisSpec axisZvtx{100, -15., 15., "z_{vtx} (cm)"};
    const AxisSpec axisMult{100, 0., 180., "N_{ch}"};
    const AxisSpec axisEta{30, -1.5, 1.5, "#eta"};
    const AxisSpec axisPhi{60, 0., 2. * TMath::Pi(), "#phi"};
    const AxisSpec axisDphi{30, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(), "#Delta#phi"};
    const AxisSpec axisDeta{32, -1.6, 1.6, "#Delta#eta"};
    // Reserved for a future rapidity-differentiated analysis -- unused
    // while the PID correlations run pseudorapidity-only (see
    // eventshape-pid-pseudorapidity-only auto-memory, 2026-09-15).
    [[maybe_unused]] const AxisSpec axisDy{32, -1.6, 1.6, "#Delta y"};
    const AxisSpec axisPt{nBinsPt, 0., 10., "p_{T} (GeV/c)"};
    const AxisSpec axisSph{100, 0., 1., "S_{T}"};
    const AxisSpec axisPt2{nBinsPt2, 0., 5., "p_{T} (GeV/c)"};

    // ----------------------------------------------------------------
    //  Lambda comuna: creeaza un set complet de histograme pentru un
    //  "canal" (MC sau Reco), cu toate numele suffixate corespunzator.
    //  Primeste referinte la toti membrii specifici canalului.
    // ----------------------------------------------------------------
    auto bookChannel = [&](const char* suf,
                           std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses>& hSE_,
                           std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses>& hME_,
                           std::array<std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses>, nPidCorr>& hSEpid_,
                           std::array<std::array<std::array<std::shared_ptr<TH2>, nMultBins>, nSphClasses>, nPidCorr>& hMEpid_,
                           std::shared_ptr<TH2>& hNtrig_,
                           std::array<std::shared_ptr<TH2>, nPidCorr>& hNtrigPID_,
                           std::array<std::shared_ptr<TH1>, nMultBins>& hSphMult_) {
      // QA
      // evSel bins: 1=read 2=zvtx 3=multBin 4=zvtxBin 5=leadPt 6=ST_valid 7=sphClass 8=SE_filled
      auto hQA = histos.add<TH1>(Form("evSel_%s", suf), Form("Event selection (%s)", suf), HistType::kTH1D, {{8, 0.5, 8.5}});
      hQA->GetXaxis()->SetBinLabel(1, "Events read");
      hQA->GetXaxis()->SetBinLabel(2, "|z_{vtx}|<10");
      hQA->GetXaxis()->SetBinLabel(3, "multBin OK");
      hQA->GetXaxis()->SetBinLabel(4, "zvtxBin OK");
      hQA->GetXaxis()->SetBinLabel(5, "leadPt #in [min,max]");
      hQA->GetXaxis()->SetBinLabel(6, "S_{T} valid");
      hQA->GetXaxis()->SetBinLabel(7, "sphClass OK");
      hQA->GetXaxis()->SetBinLabel(8, "SE filled");

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
      for (int im = 0; im < nMultBins; ++im) {
        hSphMult_[im] = histos.add<TH1>(
          Form("sphericity_mult%d_%s", im, suf),
          Form("S_T [%s] (%s);S_{T};events", multLabel(im).c_str(), suf),
          HistType::kTH1D, {axisSph});
      }

      // Trigger counter TH2: x = multBin, y = sphClass
      hNtrig_ = histos.add<TH2>(
        Form("hNtrig_%s", suf),
        Form("N_{trig} (%s);mult bin;sph class", suf),
        HistType::kTH2D,
        {{nMultBins, -0.5, static_cast<double>(nMultBins) - 0.5},
         {nSphClasses, -0.5, static_cast<double>(nSphClasses) - 0.5}});
      for (int im = 0; im < nMultBins; ++im) {
        hNtrig_->GetXaxis()->SetBinLabel(im + 1, multLabel(im).c_str());
      }
      for (int ic = 0; ic < nSphClasses; ++ic) {
        hNtrig_->GetYaxis()->SetBinLabel(ic + 1, kSphLabels[ic]);
      }

      // SE/ME TH2D: [sphClass][multBin]
      for (int ic = 0; ic < nSphClasses; ++ic) {
        for (int im = 0; im < nMultBins; ++im) {
          hSE_[ic][im] = histos.add<TH2>(
            Form("hSE_%s_sph%d_mult%d", suf, ic, im),
            Form("SE [%s][%s] (%s);#Delta#phi;#Delta#eta", kSphLabels[ic], multLabel(im).c_str(), suf),
            HistType::kTH2D, {axisDphi, axisDeta});
          hME_[ic][im] = histos.add<TH2>(
            Form("hME_%s_sph%d_mult%d", suf, ic, im),
            Form("ME [%s][%s] (%s);#Delta#phi;#Delta#eta", kSphLabels[ic], multLabel(im).c_str(), suf),
            HistType::kTH2D, {axisDphi, axisDeta});
        }
      }

      // ---- SE/ME like-sign per specie PID: [pidIdx][sphClass][multBin] ---
      // Pseudorapidity-only regime (confirmed 2026-09-15, see Claude Code
      // auto-memory eventshape-pid-pseudorapidity-only): PID correlations
      // use Delta-eta, same as the inclusive hSE_/hME_ histograms, NOT
      // Delta-y. The rapidity infrastructure (computeRapidity(), the
      // per-track y field, axisDy below) is kept for a future
      // higher-statistics rapidity-differentiated analysis -- it is
      // deliberately unused here for now, not dead code to remove.
      for (int ip = 0; ip < nPidCorr; ++ip) {
        for (int ic = 0; ic < nSphClasses; ++ic) {
          for (int im = 0; im < nMultBins; ++im) {
            hSEpid_[ip][ic][im] = histos.add<TH2>(
              Form("hSE_%s_pid%d_sph%d_mult%d", suf, ip, ic, im),
              Form("SE %s--%s [%s][%s] (%s);#Delta#varphi;#Delta#eta",
                   kPidCorrNames[ip], kPidCorrNames[ip], kSphLabels[ic], multLabel(im).c_str(), suf),
              HistType::kTH2D, {axisDphi, axisDeta});
            hMEpid_[ip][ic][im] = histos.add<TH2>(
              Form("hME_%s_pid%d_sph%d_mult%d", suf, ip, ic, im),
              Form("ME %s--%s [%s][%s] (%s);#Delta#varphi;#Delta#eta",
                   kPidCorrNames[ip], kPidCorrNames[ip], kSphLabels[ic], multLabel(im).c_str(), suf),
              HistType::kTH2D, {axisDphi, axisDeta});
          }
        }
        // Trigger counter per specie PID
        hNtrigPID_[ip] = histos.add<TH2>(
          Form("hNtrig_%s_pid%d", suf, ip),
          Form("N_{trig} %s (%s);mult bin;sph class", kPidCorrNames[ip], suf),
          HistType::kTH2D,
          {{nMultBins, -0.5, static_cast<double>(nMultBins) - 0.5},
           {nSphClasses, -0.5, static_cast<double>(nSphClasses) - 0.5}});
        for (int im = 0; im < nMultBins; ++im) {
          hNtrigPID_[ip]->GetXaxis()->SetBinLabel(im + 1, multLabel(im).c_str());
        }
        for (int ic = 0; ic < nSphClasses; ++ic) {
          hNtrigPID_[ip]->GetYaxis()->SetBinLabel(ic + 1, kSphLabels[ic]);
        }
      }
    };

    // Book histograms for MC (generated) channel
    bookChannel("MC",
                hSE_MC, hME_MC, hSEpid_MC, hMEpid_MC,
                hNtrig_MC, hNtrigPID_MC, hSphMult_MC);

    // Book histograms for Reco (reconstructed) channel
    bookChannel("Reco",
                hSE_Reco, hME_Reco, hSEpid_Reco, hMEpid_Reco,
                hNtrig_Reco, hNtrigPID_Reco, hSphMult_Reco);

    // [DEBUG, TEMPORARY] track-cut diagnostic counter for processReco
    hTrackCutDebug_Reco = histos.add<TH1>(
      "hTrackCutDebug_Reco",
      "Track cut flow (Reco, DEBUG);cut;tracks",
      HistType::kTH1D, {{6, 0.5, 6.5}});
    hTrackCutDebug_Reco->GetXaxis()->SetBinLabel(1, "seen");
    hTrackCutDebug_Reco->GetXaxis()->SetBinLabel(2, "passed isGlobalTrackWoDCA");
    hTrackCutDebug_Reco->GetXaxis()->SetBinLabel(3, "passed DCAxy");
    hTrackCutDebug_Reco->GetXaxis()->SetBinLabel(4, "passed DCAz");
    hTrackCutDebug_Reco->GetXaxis()->SetBinLabel(5, "passed eta");
    hTrackCutDebug_Reco->GetXaxis()->SetBinLabel(6, "passed pt (selected)");

    // PID confusion matrix (rapid contamination estimate, 2026-08-14) --
    // axes: pT, true species (PDG), tagged species (pidFromNSigma).
    // Species axis order matches PidSpecies: 0=unid,1=pion,2=kaon,3=proton.
    const AxisSpec axisPidSpecies{static_cast<int>(kNPidSpecies), -0.5, static_cast<double>(kNPidSpecies) - 0.5, "species"};
    hPidConfusion_MC = histos.add<TH3>(
      "hPidConfusion_MC",
      "True vs. tagged PID species (MC);p_{T} (GeV/c);true species;tagged species",
      HistType::kTH3D, {axisPt2, axisPidSpecies, axisPidSpecies});
    for (int ip = 0; ip < static_cast<int>(kNPidSpecies); ++ip) {
      hPidConfusion_MC->GetYaxis()->SetBinLabel(ip + 1, kPidNames[ip]);
      hPidConfusion_MC->GetZaxis()->SetBinLabel(ip + 1, kPidNames[ip]);
    }
  }

  // ----------------------------------------------------------------
  //  processMC
  // ----------------------------------------------------------------
  void processMC(aod::McCollision const& mcCollision,
                 aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("evSel_MC"), 1);

    const double zvtx = mcCollision.posZ();
    if (std::fabs(zvtx) > maxZvtx) {
      return;
    }
    histos.fill(HIST("evSel_MC"), 2);
    histos.fill(HIST("Zvertex_MC"), zvtx);

    // ---- one-pass: build selTracks + find leading ---------------
    std::vector<TrackSimple> selTracks;
    selTracks.reserve(kMultEdges.back() - 1);
    int leadIdx = -1;
    double pTlead = -1., phiLead = 0., etaLead = 0.;

    for (const auto& p : mcParticles) {
      if (!p.isPhysicalPrimary()) {
        continue;
      }
      auto* pdgPtr = pdg->GetParticle(p.pdgCode());
      if (!pdgPtr || std::fabs(pdgPtr->Charge()) < 0.1) {
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
      // kPidUnidentified via pidFromPDG) when the rapidity cut fails, so the
      // track contributes to hSE/hME but not hSEpid/hMEpid.
      int pdgForPID = p.pdgCode();
      if (applyRapidityCutPID && pidFromPDG(pdgForPID) != kPidUnidentified && std::fabs(rap) > maxRapPID) {
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

    const int mBin = o2::analysis::findBin(&kMultEdges, nch);
    if (mBin < 0) {
      return;
    }
    histos.fill(HIST("evSel_MC"), 3);

    const int zBin = o2::analysis::findBin(&kZvtxEdges, zvtx);
    if (zBin < 0) {
      return;
    }
    histos.fill(HIST("evSel_MC"), 4);

    // ---- leading particle pT cut [minLeadPt, maxLeadPt] -----------
    if (leadIdx < 0) {
      return;
    }
    if (pTlead < minLeadPt || pTlead > maxLeadPt) {
      return;
    }
    histos.fill(HIST("evSel_MC"), 5);
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
    histos.fill(HIST("evSel_MC"), 6);
    histos.fill(HIST("sphericity_MC"), sph);
    hSphMult_MC[mBin]->Fill(sph);

    const int sphCls = sphClass(sph); // 0,1,2 or -1
    if (sphCls < 0) {
      return;
    }
    histos.fill(HIST("evSel_MC"), 7);

    const PidSpecies pidLead = pidFromPDG(selTracks[leadIdx].pdg);

    // ---- trigger counter ----------------------------------------
    hNtrig_MC->Fill(static_cast<double>(mBin), static_cast<double>(sphCls));
    hNtrig_MC->Fill(static_cast<double>(mBin), static_cast<double>((nSphClasses - 1))); // "all"

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
        continue; // asociat pT < leading pT
      }
      if (assoc.pt < minAssocCorr || assoc.pt > maxAssocCorr) {
        continue;
      }

      histos.fill(HIST("ptAssocHistogram_MC"), assoc.pt);

      const PidSpecies pidAssoc = pidFromPDG(assoc.pdg);

      const double dphi = RecoDecay::constrainAngle(phiLead - assoc.phi, -o2::constants::math::PIHalf);
      const double deta = etaLead - assoc.eta;

      // TH2D SE: specific sph class (uses Delta-eta)
      hSE_MC[sphCls][mBin]->Fill(dphi, deta);
      // TH2D SE: "all" class (uses Delta-eta)
      hSE_MC[nSphClasses - 1][mBin]->Fill(dphi, deta);

      // SE like-sign: fill solo se trigger e associato sono la stessa specie
      // Uses Delta-eta (pseudorapidity-only regime, see
      // eventshape-pid-pseudorapidity-only auto-memory)
      int pidLeadIdx = pidSpeciesToIdx(pidLead);
      int pidAssocIdx = pidSpeciesToIdx(pidAssoc);
      if (pidLeadIdx >= 0 && pidLeadIdx == pidAssocIdx) {
        hSEpid_MC[pidLeadIdx][sphCls][mBin]->Fill(dphi, deta);
        hSEpid_MC[pidLeadIdx][nSphClasses - 1][mBin]->Fill(dphi, deta);
      }

      seHasPairs = true;
    }
    if (seHasPairs) {
      histos.fill(HIST("evSel_MC"), 8);
    }

    // ---- hNtrig_MC per specie PID leading --------------------------
    int pidLeadIdx = pidSpeciesToIdx(pidLead);
    if (pidLeadIdx >= 0) {
      hNtrigPID_MC[pidLeadIdx]->Fill(static_cast<double>(mBin), static_cast<double>(sphCls));
      hNtrigPID_MC[pidLeadIdx]->Fill(static_cast<double>(mBin), static_cast<double>((nSphClasses - 1)));
    }

    // ================================================================
    //  Update pools AFTER SE
    // ================================================================
    // Pool "all" — specific sph
    {
      auto& pool = eventPools_MC[sphCls][zBin][mBin];
      pool.push_back(selTracks);
      if (static_cast<int>(pool.size()) > mixPoolSize) {
        pool.pop_front();
      }
    }
    // Pool "all" — all sph
    {
      auto& poolAll = eventPools_MC[nSphClasses - 1][zBin][mBin];
      poolAll.push_back(selTracks);
      if (static_cast<int>(poolAll.size()) > mixPoolSize) {
        poolAll.pop_front();
      }
    }

    // Pool per specie PID — contiene solo i track di quella specie
    for (int ip = 0; ip < nPidCorr; ++ip) {
      // Costruisci un vettore con solo i track della specie ip
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
        auto& pool = eventPoolsPID_MC[ip][sphCls][zBin][mBin];
        pool.push_back(pidTracks);
        if (static_cast<int>(pool.size()) > mixPoolSize) {
          pool.pop_front();
        }
      }
      // "all" sph pool
      {
        auto& pool = eventPoolsPID_MC[ip][nSphClasses - 1][zBin][mBin];
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
      const auto& pool = eventPools_MC[poolIdx][zBin][mBin];
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
          hME_MC[poolIdx][mBin]->Fill(dphi, deta);
        }
      }
    };

    doMixing(sphCls);
    doMixing(nSphClasses - 1);

    // ================================================================
    //  MIXED-EVENT like-sign per specie PID
    //  Trigger = leading della specie pidLeadIdx
    //  Associati = track della stessa specie dal pool PID
    //  Uses Delta-eta (pseudorapidity-only regime, see
    //  eventshape-pid-pseudorapidity-only auto-memory)
    // ================================================================
    if (pidLeadIdx >= 0) {
      auto doMixingPID = [&](int poolIdx) {
        const auto& pool = eventPoolsPID_MC[pidLeadIdx][poolIdx][zBin][mBin];
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
            hMEpid_MC[pidLeadIdx][poolIdx][mBin]->Fill(dphi, detaPID);
          }
        }
      };
      doMixingPID(sphCls);
      doMixingPID(nSphClasses - 1);
    }
  }
  PROCESS_SWITCH(myExampleTask, processMC, "Process MC truth", true);

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
    histos.fill(HIST("evSel_Reco"), 1);

    // ---- event selection ----
    if (!col.sel8()) {
      return;
    }
    const double zvtx = col.posZ();
    if (std::fabs(zvtx) > maxZvtx) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), 2);
    histos.fill(HIST("Zvertex_Reco"), zvtx);

    // ---- one-pass: build selTracks + find leading ---------------
    std::vector<TrackSimple> selTracks;
    selTracks.reserve(kMultEdges.back() - 1);
    int leadIdx = -1;
    double pTlead = -1., phiLead = 0., etaLead = 0.;

    for (const auto& track : tracks) {
      hTrackCutDebug_Reco->Fill(1); // seen

      // Quality selection: global track without DCA (apply DCA separately)
      if (!track.isGlobalTrackWoDCA()) {
        continue;
      }
      hTrackCutDebug_Reco->Fill(2); // passed isGlobalTrackWoDCA

      // DCA cuts
      if (std::fabs(track.dcaXY()) > maxDcaXY) {
        continue;
      }
      hTrackCutDebug_Reco->Fill(3); // passed DCAxy
      if (std::fabs(track.dcaZ()) > maxDcaZ) {
        continue;
      }
      hTrackCutDebug_Reco->Fill(4); // passed DCAz

      // Kinematic cuts (same as MC)
      if (std::fabs(track.eta()) > maxEta) {
        continue;
      }
      hTrackCutDebug_Reco->Fill(5); // passed eta
      if (track.pt() < minAssocPt) {
        continue;
      }
      hTrackCutDebug_Reco->Fill(6); // passed pt (selected)

      // PID from nSigma (combined TPC+TOF, ALICE 2511.10399, unless the
      // legacy flags below switch to pre-2026-08-11 gating/selection)
      PidSpecies pid = pidFromNSigma(track, maxNSigmaPID, maxNSigmaPIDAmbig,
                                     maxNSigmaTPC, maxNSigmaTOF,
                                     useLegacyNSigmaGating, useLegacyPIDSelection);
      double mass = pidMass(pid);
      double rap = computeRapidity(track.pt(), track.eta(), mass);

      // |y| < maxRapPID cut applies ONLY to PID-identified species -- see
      // matching comment in processMC. Track stays in selTracks either way
      // (all-species channel unaffected); only the PID assignment is reset.
      if (applyRapidityCutPID && pid != kPidUnidentified && std::fabs(rap) > maxRapPID) {
        pid = kPidUnidentified;
      }

      // pdg = 0 for reco (no MC truth); store sign for future like/unlike-sign
      int idx = static_cast<int>(selTracks.size());
      selTracks.push_back({static_cast<float>(track.phi()), static_cast<float>(track.eta()), static_cast<float>(track.pt()), 0});
      // Store PID info in pdg field as signed pidIdx for later use:
      //   sign(track.sign()) * (211 for pion, 321 for kaon, 2212 for proton, 0 for unid)
      int pdgLike = 0;
      if (pid == kPidPion) {
        pdgLike = 211;
      }
      if (pid == kPidKaon) {
        pdgLike = 321;
      }
      if (pid == kPidProton) {
        pdgLike = 2212;
      }
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

    const int mBin = o2::analysis::findBin(&kMultEdges, nch);
    if (mBin < 0) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), 3);

    const int zBin = o2::analysis::findBin(&kZvtxEdges, zvtx);
    if (zBin < 0) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), 4);

    // ---- leading particle pT cut [minLeadPt, maxLeadPt] -----------
    if (leadIdx < 0) {
      return;
    }
    if (pTlead < minLeadPt || pTlead > maxLeadPt) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), 5);
    histos.fill(HIST("ptLeadHistogram_Reco"), pTlead);

    // ---- sphericity ---------------------------------------------
    const double sph = computeSphericity(selTracks);
    if (sph >= 0.) {
      histos.fill(HIST("sphericity_beforeCuts_Reco"), sph);
    }
    if (sph < 0.) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), 6);
    histos.fill(HIST("sphericity_Reco"), sph);
    hSphMult_Reco[mBin]->Fill(sph);

    const int sphCls = sphClass(sph);
    if (sphCls < 0) {
      return;
    }
    histos.fill(HIST("evSel_Reco"), 7);

    const PidSpecies pidLead = pidFromPDG(selTracks[leadIdx].pdg);

    // ---- trigger counter ----------------------------------------
    hNtrig_Reco->Fill(static_cast<double>(mBin), static_cast<double>(sphCls));
    hNtrig_Reco->Fill(static_cast<double>(mBin), static_cast<double>((nSphClasses - 1)));

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
      hSE_Reco[sphCls][mBin]->Fill(dphi, deta);
      hSE_Reco[nSphClasses - 1][mBin]->Fill(dphi, deta);

      // PID SE — same species, Delta-eta (pseudorapidity-only regime, see
      // eventshape-pid-pseudorapidity-only auto-memory)
      int pidLeadIdx = pidSpeciesToIdx(pidLead);
      int pidAssocIdx = pidSpeciesToIdx(pidAssoc);
      if (pidLeadIdx >= 0 && pidLeadIdx == pidAssocIdx) {
        hSEpid_Reco[pidLeadIdx][sphCls][mBin]->Fill(dphi, deta);
        hSEpid_Reco[pidLeadIdx][nSphClasses - 1][mBin]->Fill(dphi, deta);
      }

      seHasPairs = true;
    }
    if (seHasPairs) {
      histos.fill(HIST("evSel_Reco"), 8);
    }

    // ---- hNtrig_Reco per specie PID leading --------------------------
    int pidLeadIdx = pidSpeciesToIdx(pidLead);
    if (pidLeadIdx >= 0) {
      hNtrigPID_Reco[pidLeadIdx]->Fill(static_cast<double>(mBin), static_cast<double>(sphCls));
      hNtrigPID_Reco[pidLeadIdx]->Fill(static_cast<double>(mBin), static_cast<double>((nSphClasses - 1)));
    }

    // ================================================================
    //  Update pools AFTER SE (identical to processMC)
    // ================================================================
    {
      auto& pool = eventPools_Reco[sphCls][zBin][mBin];
      pool.push_back(selTracks);
      if (static_cast<int>(pool.size()) > mixPoolSize) {
        pool.pop_front();
      }
    }
    {
      auto& poolAll = eventPools_Reco[nSphClasses - 1][zBin][mBin];
      poolAll.push_back(selTracks);
      if (static_cast<int>(poolAll.size()) > mixPoolSize) {
        poolAll.pop_front();
      }
    }

    for (int ip = 0; ip < nPidCorr; ++ip) {
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
        auto& pool = eventPoolsPID_Reco[ip][sphCls][zBin][mBin];
        pool.push_back(pidTracks);
        if (static_cast<int>(pool.size()) > mixPoolSize) {
          pool.pop_front();
        }
      }
      {
        auto& pool = eventPoolsPID_Reco[ip][nSphClasses - 1][zBin][mBin];
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
      const auto& pool = eventPools_Reco[poolIdx][zBin][mBin];
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
          hME_Reco[poolIdx][mBin]->Fill(dphi, deta);
        }
      }
    };

    doMixing(sphCls);
    doMixing(nSphClasses - 1);

    // ================================================================
    //  MIXED-EVENT PID (Delta-eta, pseudorapidity-only regime, see
    //  eventshape-pid-pseudorapidity-only auto-memory)
    // ================================================================
    if (pidLeadIdx >= 0) {
      auto doMixingPID = [&](int poolIdx) {
        const auto& pool = eventPoolsPID_Reco[pidLeadIdx][poolIdx][zBin][mBin];
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
            hMEpid_Reco[pidLeadIdx][poolIdx][mBin]->Fill(dphi, detaPID);
          }
        }
      };
      doMixingPID(sphCls);
      doMixingPID(nSphClasses - 1);
    }
  }
  PROCESS_SWITCH(myExampleTask, processReco, "Process reconstructed data", false);

  // ----------------------------------------------------------------
  //  processMCConfusion -- "rapid" PID contamination estimate, requested
  //  2026-08-14. Fills hPidConfusion_MC (pT, true species, tagged
  //  species) on MC only -- needs BOTH the truth PDG (via McTrackLabels)
  //  AND the reconstructed track's nSigma values on the SAME track,
  //  which neither processMC (truth only) nor processReco (no truth
  //  link) has access to, hence a dedicated process function. Reports
  //  contamination as a systematic; does NOT correct anything (that's
  //  the separate, not-yet-started "riguros" unfolding approach -- see
  //  STATUS.md). Off by default, like the other MC-only diagnostics.
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
      if (!pdgPtr || std::fabs(pdgPtr->Charge()) < 0.1) {
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
      // effect on this matrix (found 2026-08-14, comparing Phase 1
      // rapidityOnly against baseline: identical confusion-matrix numbers
      // despite different processReco yields, before this fix).
      if (applyRapidityCutPID && taggedPid != kPidUnidentified) {
        double tagMass = pidMass(taggedPid);
        double tagRap = computeRapidity(track.pt(), track.eta(), tagMass);
        if (std::fabs(tagRap) > maxRapPID) {
          taggedPid = kPidUnidentified;
        }
      }

      histos.fill(HIST("hPidConfusion_MC"), track.pt(), static_cast<double>(truePid), static_cast<double>(taggedPid));
    }
  }
  PROCESS_SWITCH(myExampleTask, processMCConfusion, "Rapid PID contamination estimate: true (PDG) vs tagged (nSigma) species confusion matrix, MC only", false);

}; // end struct

// ----------------------------------------------------------------
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<myExampleTask>(cfgc)};
}
