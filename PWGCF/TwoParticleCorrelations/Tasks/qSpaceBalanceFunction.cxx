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

/// \file qSpaceBalanceFunction.cxx
/// \brief K+K- pair densities in Q_inv, Q_long, Q_out and Q_side, the inputs of q-space balance functions
/// \author Lakshitha Fernando <lakshithafernando@wayne.edu>
///
/// Reads the tables of dpt-dpt-filter, which must run in the same workflow with only kaons enabled,
/// so that trackacceptedid is 0 for K+ and 1 for K-.
/// Pairs are ORDERED, as in dpt-dpt-correlations: every particle takes the role of particle 1 in turn.
/// Each of ++, +-, -+, -- is filled; an unordered like-sign pair therefore enters twice, an unordered
/// opposite-sign pair once in +- and once in -+. OS = (+-) + (-+) and LS = (++) + (--) are comparable.
/// Every pair is also classified by the MC mothers of the two kaons: phi pair (same mother, a phi(1020)),
/// other common mother, or no common mother. At detector level (processRecLevelMcNotStored) the mothers come from the
/// MC particles matched to the tracks, the true species gives the pair purity, and the true momenta give the
/// resolution Q_rec - Q_true. Detector-level histograms are stored under rec/, generator-level ones at the top level.
///
/// Balance-function inputs (folder bf/, rec/bf/): events are sorted into classes of (z_vtx, multiplicity percentile)
/// with the binning policy of the event mixing, class = iz + nz * imult. Per class the task stores
///  - same-event ordered pairs i != j (both orders, as above) vs Q: summed over events, their integral is
///    sum N_a (N_b - delta_ab), i.e. factorial moments for like-sign pairs;
///  - mixed-event ordered pairs: each unordered pair of events (e1, e2) of one class is filled in BOTH orders,
///    (a from e1, b from e2) and (a from e2, b from e1), so the reference per event is mixed / (2 x number of event pairs),
///    whose integral is <N_a><N_b>, the uncorrelated product rho1 x rho1;
///  - events, event pairs, kaons per charge and sum N_a (N_b - delta_ab), for the normalisation and the integral checks.
/// Then rho2 = same / N_events, rho1rho1 = mixed / (2 N_eventPairs), C2 = rho2 - rho1rho1, A2^ab = C2^ab / <N_b>, done offline.

#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TwoParticleCorrelations/Core/QSpaceKinematics.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
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
#include <TH2.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
constexpr int NCharges = 2;
const std::array<std::string, NCharges> chargeTags{"P", "M"};
enum IdCheckBin {
  kIdCheckGood = 1,
  kIdCheckNotKaonId,
  kIdCheckWrongPdg,
  kIdCheckNoMcTruth
};
/// Pair origin classes from the MC mother of the two kaons
enum PairOrigin {
  kPairPhi = 0,
  kPairOtherCommonMother,
  kPairNoCommonMother,
  NPairOrigins
};
const std::array<std::string, NPairOrigins> pairOriginTags{"phi", "otherMother", "noCommonMother"};
/// MC truth of a detector-level pair
enum PairTruthBin {
  kPairBothTrueKaons = 1,
  kPairMisidentified,
  kPairUnmatched,
  kPairSameMcParticle
};
/// PDG codes above this value are hadrons; below are quarks, leptons, gauge bosons and generator-internal codes
static constexpr int KPdgHadronThreshold = 100; // o2-linter: disable=pdg/explicit-code (not a PDG code)
/// Generator status codes of beam and event-setup entries (PYTHIA 11-19) and the HepMC beam status (4).
/// Same criterion as PWGCF/TwoParticleCorrelations/Tasks/particleOriginAnalysis.cxx (V. Gonzalez)
static constexpr int KBeamStatusCodeMin = 11;
static constexpr int KBeamStatusCodeMax = 19;
static constexpr int KHepMCBeamStatusCode = 4;
/// Mother of a particle: global index (-1 if none) and PDG code (0 if none)
struct MotherInfo {
  int64_t index = -1;
  int pdgCode = 0;
};
/// MC truth of one accepted particle (generator level) or track (detector level), computed once per collision
struct TruthInfo {
  bool hasMc = false;   ///< generator level: always; detector level: the track has an MC label
  int64_t mcIndex = -1; ///< global index of the MC particle
  int pdgCode = 0;      ///< true PDG code
  float pt = 0.f;       ///< true kinematics
  float eta = 0.f;
  float phi = 0.f;
  MotherInfo mother;
};
/// Combination index of the per-class factorial moments sum N_a (N_b - delta_ab): PP, PM, MP, MM
constexpr int NChargeCombinations = NCharges * NCharges;
/// Balance-function inputs of one level, per event class (z_vtx, multiplicity)
struct BfHistos {
  using Histo2DPerCharge = std::array<std::shared_ptr<TH2>, NCharges>;
  using Histo2DPerPair = std::array<Histo2DPerCharge, NCharges>;
  std::shared_ptr<TH1> fhEventsPerClass;          ///< accepted events per class
  std::shared_ptr<TH1> fhEventsOutsideClasses;    ///< events outside the class binning (not used in bf/)
  std::shared_ptr<TH1> fhMixedEventPairsPerClass; ///< unordered pairs of events mixed, per class
  std::shared_ptr<TH2> fhKaonsPerClass;           ///< sum over events of N_a, per class and charge
  std::shared_ptr<TH2> fhFactorialPerClass;       ///< sum over events of N_a (N_b - delta_ab), per class and charge combination
  Histo2DPerPair fhSameQinv, fhSameQlong, fhSameQout, fhSameQside;
  Histo2DPerPair fhMixedQinv, fhMixedQlong, fhMixedQout, fhMixedQside;
};
/// All histograms of one level (generator or detector)
struct PairHistos {
  using HistoPerCharge = std::array<std::shared_ptr<TH1>, NCharges>;
  using Histo2DPerCharge = std::array<std::shared_ptr<TH2>, NCharges>;
  using Histo2DPerPair = std::array<Histo2DPerCharge, NCharges>;
  using HistoPerPair = std::array<HistoPerCharge, NCharges>;
  bool booked = false;
  bool withResolution = false;
  std::shared_ptr<TH2> fhEvents;
  std::shared_ptr<TH1> fhIdCheck;
  Histo2DPerCharge fhPtVsCent;
  HistoPerCharge fhEta, fhPhi, fhPtFromPhi;
  Histo2DPerPair fhQinv, fhQlong, fhQout, fhQside;
  HistoPerPair fhMass, fhSumRuleResidual, fhMassRuleResidual;
  /// pair histograms split by pair origin: [origin][charge1][charge2]
  std::array<Histo2DPerPair, NPairOrigins> fhQinvOrigin, fhQlongOrigin, fhQoutOrigin, fhQsideOrigin;
  std::shared_ptr<TH1> fhCommonMotherPdg;
  /// detector level only: pair truth status and resolution (Q_true, Q_rec - Q_true) for pairs of two true kaons
  HistoPerPair fhPairTruth;
  Histo2DPerPair fhResQinv, fhResQlong, fhResQout, fhResQside;
  BfHistos bf;
};
} // namespace

struct QSpaceBalanceFunction {
  ConfigurableAxis cfgBinsQinv{"cfgBinsQinv", {2500, 0., 5.}, "Q_inv binning (GeV/c)"};
  ConfigurableAxis cfgBinsQComponent{"cfgBinsQComponent", {1000, -5., 5.}, "Q_long, Q_out, Q_side binning (GeV/c)"};
  ConfigurableAxis cfgBinsMass{"cfgBinsMass", {6000, 0., 6.}, "Pair invariant mass binning (GeV/c^2)"};
  ConfigurableAxis cfgBinsPt{"cfgBinsPt", {18, 0.2, 2.}, "pT binning (GeV/c)"};
  ConfigurableAxis cfgBinsEta{"cfgBinsEta", {16, -0.8, 0.8}, "Pseudorapidity binning"};
  ConfigurableAxis cfgBinsResidual{"cfgBinsResidual", {200, -1.e-10, 1.e-10}, "Residual binning (GeV^2/c^2)"};
  ConfigurableAxis cfgBinsQinvTrue{"cfgBinsQinvTrue", {50, 0., 1.}, "true Q_inv binning of the resolution maps (GeV/c)"};
  ConfigurableAxis cfgBinsQComponentTrue{"cfgBinsQComponentTrue", {50, -1., 1.}, "true Q component binning of the resolution maps (GeV/c)"};
  ConfigurableAxis cfgBinsQResolution{"cfgBinsQResolution", {200, -0.05, 0.05}, "Q_rec - Q_true binning (GeV/c)"};
  ConfigurableAxis cfgMixBinsZ{"cfgMixBinsZ", {VARIABLE_WIDTH, -7., -5., -3., -1., 1., 3., 5., 7.}, "z_vtx edges of the event classes (cm), same and mixed events"};
  ConfigurableAxis cfgMixBinsMult{"cfgMixBinsMult", {VARIABLE_WIDTH, 0., 1., 2., 5., 10., 15., 20., 30., 40., 50., 60., 70., 80., 90., 100.1}, "multiplicity percentile edges of the event classes, same and mixed events (the filter accepts 0-100%)"};
  Configurable<int> cfgMixDepth{"cfgMixDepth", 10, "each event is mixed with the next cfgMixDepth events of its class"};

  using BinningZMultGen = ColumnBinningPolicy<aod::mccollision::PosZ, aod::dptdptfilter::DptDptCFCollisionCentMult>;
  using BinningZMultRec = ColumnBinningPolicy<aod::collision::PosZ, aod::dptdptfilter::DptDptCFCollisionCentMult>;
  /// Plain members (not std::optional: the framework counts task members by brace initialisation, which std::optional breaks).
  /// Built here with the configurable defaults and rebuilt in init from the configured values; under/overflow give class -1.
  BinningZMultGen binningGen{{std::vector<double>{VARIABLE_WIDTH, -7., 7.}, std::vector<double>{VARIABLE_WIDTH, 0., 100.1}}, true};
  BinningZMultRec binningRec{{std::vector<double>{VARIABLE_WIDTH, -7., 7.}, std::vector<double>{VARIABLE_WIDTH, 0., 100.1}}, true};
  SliceCache cache;

  using GenCollisions = soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>;
  using GenParticles = soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>>;
  using RecCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>;
  using RecTracks = soa::Filtered<soa::Join<aod::Tracks, aod::DptDptCFTracksInfo>>;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  PairHistos genHistos; ///< generator level, top-level folders (as before the detector level existed)
  PairHistos recHistos; ///< detector level, folder rec/

  Filter onlyAcceptedCollisions = (aod::dptdptfilter::collisionaccepted == uint8_t(true));
  Filter onlyAcceptedTracks = (aod::dptdptfilter::trackacceptedid >= int8_t(0));

  /// Books the histograms of one level under the folder prefix ("" or "rec/")
  void bookHistos(PairHistos& histos, const std::string& prefix, bool withResolution)
  {
    const AxisSpec axisCent{std::vector<double>{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.1}, "multiplicity percentile"};
    const AxisSpec axisZ{28, -7., 7., "z_{vtx} (cm)"};
    const AxisSpec axisPt{cfgBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisEta{cfgBinsEta, "#eta"};
    const AxisSpec axisPhi{72, 0., constants::math::TwoPI, "#varphi (rad)"};
    const AxisSpec axisQinv{cfgBinsQinv, "Q_{inv} (GeV/#it{c})"};
    const AxisSpec axisQlong{cfgBinsQComponent, "Q_{long} (GeV/#it{c})"};
    const AxisSpec axisQout{cfgBinsQComponent, "Q_{out} (GeV/#it{c})"};
    const AxisSpec axisQside{cfgBinsQComponent, "Q_{side} (GeV/#it{c})"};
    const AxisSpec axisMass{cfgBinsMass, "M_{KK} (GeV/#it{c}^{2})"};
    const AxisSpec axisResidual{cfgBinsResidual, "residual (GeV^{2}/#it{c}^{2})"};

    histos.withResolution = withResolution;
    histos.fhEvents = registry.add<TH2>((prefix + "events/hZVsCent").c_str(), "analysed events", kTH2D, {axisCent, axisZ});
    histos.fhIdCheck = registry.add<TH1>((prefix + "singles/hIdCheck").c_str(), "particle id vs MC truth", kTH1D, {{4, 0.5, 4.5}});
    histos.fhIdCheck->GetXaxis()->SetBinLabel(kIdCheckGood, "kaon, id agrees with true PDG");
    histos.fhIdCheck->GetXaxis()->SetBinLabel(kIdCheckNotKaonId, "id is not 0 or 1");
    histos.fhIdCheck->GetXaxis()->SetBinLabel(kIdCheckWrongPdg, "id disagrees with true PDG");
    histos.fhIdCheck->GetXaxis()->SetBinLabel(kIdCheckNoMcTruth, "no MC label");
    for (int charge1 = 0; charge1 < NCharges; ++charge1) {
      const std::string& tag1 = chargeTags[charge1];
      histos.fhPtVsCent[charge1] = registry.add<TH2>((prefix + "singles/hPtVsCent_" + tag1).c_str(), "kaon #it{p}_{T}", kTH2D, {axisCent, axisPt});
      histos.fhEta[charge1] = registry.add<TH1>((prefix + "singles/hEta_" + tag1).c_str(), "kaon #eta", kTH1D, {axisEta});
      histos.fhPhi[charge1] = registry.add<TH1>((prefix + "singles/hPhi_" + tag1).c_str(), "kaon #varphi", kTH1D, {axisPhi});
      histos.fhPtFromPhi[charge1] = registry.add<TH1>((prefix + "singles/hPtFromPhi_" + tag1).c_str(), "kaons whose mother is a #phi(1020)", kTH1D, {axisPt});
      for (int charge2 = 0; charge2 < NCharges; ++charge2) {
        const std::string tag = tag1 + chargeTags[charge2];
        histos.fhQinv[charge1][charge2] = registry.add<TH2>((prefix + "pairs/hQinv_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQinv});
        histos.fhQlong[charge1][charge2] = registry.add<TH2>((prefix + "pairs/hQlong_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQlong});
        histos.fhQout[charge1][charge2] = registry.add<TH2>((prefix + "pairs/hQout_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQout});
        histos.fhQside[charge1][charge2] = registry.add<TH2>((prefix + "pairs/hQside_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQside});
        histos.fhMass[charge1][charge2] = registry.add<TH1>((prefix + "pairs/hMass_" + tag).c_str(), "ordered pairs", kTH1F, {axisMass});
        histos.fhSumRuleResidual[charge1][charge2] = registry.add<TH1>((prefix + "checks/hSumRuleResidual_" + tag).c_str(), "Q_{long}^{2}+Q_{out}^{2}+Q_{side}^{2}-Q_{inv}^{2}", kTH1D, {axisResidual});
        histos.fhMassRuleResidual[charge1][charge2] = registry.add<TH1>((prefix + "checks/hMassRuleResidual_" + tag).c_str(), "M^{2}-4m^{2}-Q_{inv}^{2}", kTH1D, {axisResidual});
        for (int origin = 0; origin < NPairOrigins; ++origin) {
          const std::string originTag = tag + "_" + pairOriginTags[origin];
          histos.fhQinvOrigin[origin][charge1][charge2] = registry.add<TH2>((prefix + "pairsByOrigin/hQinv_" + originTag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQinv});
          histos.fhQlongOrigin[origin][charge1][charge2] = registry.add<TH2>((prefix + "pairsByOrigin/hQlong_" + originTag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQlong});
          histos.fhQoutOrigin[origin][charge1][charge2] = registry.add<TH2>((prefix + "pairsByOrigin/hQout_" + originTag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQout});
          histos.fhQsideOrigin[origin][charge1][charge2] = registry.add<TH2>((prefix + "pairsByOrigin/hQside_" + originTag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQside});
        }
        if (withResolution) {
          const AxisSpec axisQinvTrue{cfgBinsQinvTrue, "true Q_{inv} (GeV/#it{c})"};
          const AxisSpec axisQComponentTrue{cfgBinsQComponentTrue, "true Q component (GeV/#it{c})"};
          const AxisSpec axisResolution{cfgBinsQResolution, "Q_{rec} - Q_{true} (GeV/#it{c})"};
          histos.fhPairTruth[charge1][charge2] = registry.add<TH1>((prefix + "resolution/hPairTruth_" + tag).c_str(), "ordered pairs by MC truth", kTH1D, {{4, 0.5, 4.5}});
          histos.fhPairTruth[charge1][charge2]->GetXaxis()->SetBinLabel(kPairBothTrueKaons, "both true kaons, right charge");
          histos.fhPairTruth[charge1][charge2]->GetXaxis()->SetBinLabel(kPairMisidentified, "at least one misidentified");
          histos.fhPairTruth[charge1][charge2]->GetXaxis()->SetBinLabel(kPairUnmatched, "at least one without MC label");
          histos.fhPairTruth[charge1][charge2]->GetXaxis()->SetBinLabel(kPairSameMcParticle, "both tracks from one MC particle");
          histos.fhResQinv[charge1][charge2] = registry.add<TH2>((prefix + "resolution/hResQinv_" + tag).c_str(), "pairs of two true kaons", kTH2F, {axisQinvTrue, axisResolution});
          histos.fhResQlong[charge1][charge2] = registry.add<TH2>((prefix + "resolution/hResQlong_" + tag).c_str(), "pairs of two true kaons", kTH2F, {axisQComponentTrue, axisResolution});
          histos.fhResQout[charge1][charge2] = registry.add<TH2>((prefix + "resolution/hResQout_" + tag).c_str(), "pairs of two true kaons", kTH2F, {axisQComponentTrue, axisResolution});
          histos.fhResQside[charge1][charge2] = registry.add<TH2>((prefix + "resolution/hResQside_" + tag).c_str(), "pairs of two true kaons", kTH2F, {axisQComponentTrue, axisResolution});
        }
      }
    }
    histos.fhCommonMotherPdg = registry.add<TH1>((prefix + "pairsByOrigin/hCommonMotherPdg").c_str(), "|PDG code| of the common mother, ordered pairs (overflow = |PDG| >= 10000)", kTH1D, {{10000, 0., 10000.}});
    bookBfHistos(histos.bf, prefix + "bf/", axisQinv, axisQlong, axisQout, axisQside);
    histos.booked = true;
  }

  /// Books the balance-function inputs of one level; the class axis follows the mixing binning
  void bookBfHistos(BfHistos& bf, const std::string& prefix, const AxisSpec& axisQinv, const AxisSpec& axisQlong, const AxisSpec& axisQout, const AxisSpec& axisQside)
  {
    const int nClasses = binningGen.getAllBinsCount();
    const AxisSpec axisClass{nClasses, -0.5, nClasses - 0.5, "event class (iz + n_{z} imult)"};
    bf.fhEventsPerClass = registry.add<TH1>((prefix + "hEventsPerClass").c_str(), "accepted events", kTH1D, {axisClass});
    bf.fhEventsOutsideClasses = registry.add<TH1>((prefix + "hEventsOutsideClasses").c_str(), "events outside the class binning", kTH1D, {{1, 0.5, 1.5}});
    bf.fhMixedEventPairsPerClass = registry.add<TH1>((prefix + "hMixedEventPairsPerClass").c_str(), "unordered event pairs mixed (each filled in both orders)", kTH1D, {axisClass});
    bf.fhKaonsPerClass = registry.add<TH2>((prefix + "hKaonsPerClass").c_str(), "#sum N_{a}", kTH2D, {axisClass, {NCharges, -0.5, NCharges - 0.5, "charge (0 = K^{+}, 1 = K^{-})"}});
    bf.fhFactorialPerClass = registry.add<TH2>((prefix + "hFactorialPerClass").c_str(), "#sum N_{a}(N_{b}-#delta_{ab})", kTH2D, {axisClass, {NChargeCombinations, -0.5, NChargeCombinations - 0.5, "ab (0 = PP, 1 = PM, 2 = MP, 3 = MM)"}});
    for (int charge1 = 0; charge1 < NCharges; ++charge1) {
      for (int charge2 = 0; charge2 < NCharges; ++charge2) {
        const std::string tag = chargeTags[charge1] + chargeTags[charge2];
        bf.fhSameQinv[charge1][charge2] = registry.add<TH2>((prefix + "same/hQinv_" + tag).c_str(), "same-event ordered pairs", kTH2F, {axisClass, axisQinv});
        bf.fhSameQlong[charge1][charge2] = registry.add<TH2>((prefix + "same/hQlong_" + tag).c_str(), "same-event ordered pairs", kTH2F, {axisClass, axisQlong});
        bf.fhSameQout[charge1][charge2] = registry.add<TH2>((prefix + "same/hQout_" + tag).c_str(), "same-event ordered pairs", kTH2F, {axisClass, axisQout});
        bf.fhSameQside[charge1][charge2] = registry.add<TH2>((prefix + "same/hQside_" + tag).c_str(), "same-event ordered pairs", kTH2F, {axisClass, axisQside});
        bf.fhMixedQinv[charge1][charge2] = registry.add<TH2>((prefix + "mixed/hQinv_" + tag).c_str(), "mixed-event ordered pairs, both orders of each event pair", kTH2F, {axisClass, axisQinv});
        bf.fhMixedQlong[charge1][charge2] = registry.add<TH2>((prefix + "mixed/hQlong_" + tag).c_str(), "mixed-event ordered pairs, both orders of each event pair", kTH2F, {axisClass, axisQlong});
        bf.fhMixedQout[charge1][charge2] = registry.add<TH2>((prefix + "mixed/hQout_" + tag).c_str(), "mixed-event ordered pairs, both orders of each event pair", kTH2F, {axisClass, axisQout});
        bf.fhMixedQside[charge1][charge2] = registry.add<TH2>((prefix + "mixed/hQside_" + tag).c_str(), "mixed-event ordered pairs, both orders of each event pair", kTH2F, {axisClass, axisQside});
      }
    }
  }

  void init(InitContext const&)
  {
    binningGen = BinningZMultGen{{cfgMixBinsZ.value, cfgMixBinsMult.value}, true};
    binningRec = BinningZMultRec{{cfgMixBinsZ.value, cfgMixBinsMult.value}, true};
    if (doprocessGenLevelNotStored.value || doprocessGenLevelMixedNotStored.value) {
      bookHistos(genHistos, "", false);
    }
    if (doprocessRecLevelMcNotStored.value || doprocessRecLevelMixedNotStored.value) {
      bookHistos(recHistos, "rec/", true);
    }
  }

  /// Fills ordered pairs (a from tracksA, b from tracksB) of one event class; for the same event the pair of a
  /// particle with itself is skipped, so like-sign pairs count N(N-1)
  template <bool sameEvent, typename TTracksA, typename TTracksB>
  void fillBfPairs(BfHistos& bf, int eventClass, TTracksA const& tracksA, TTracksB const& tracksB)
  {
    const double massKaon = constants::physics::MassKaonCharged;
    for (auto const& particle1 : tracksA) {
      const int charge1 = particle1.trackacceptedid();
      if (charge1 >= NCharges) {
        continue;
      }
      for (auto const& particle2 : tracksB) {
        const int charge2 = particle2.trackacceptedid();
        if (charge2 >= NCharges) {
          continue;
        }
        if constexpr (sameEvent) {
          if (particle1.globalIndex() == particle2.globalIndex()) {
            continue;
          }
        }
        const auto pair = analysis::qspace::computePairQ(particle1.pt(), particle1.eta(), particle1.phi(), massKaon,
                                                         particle2.pt(), particle2.eta(), particle2.phi(), massKaon);
        if constexpr (sameEvent) {
          bf.fhSameQinv[charge1][charge2]->Fill(eventClass, pair.qInv);
          bf.fhSameQlong[charge1][charge2]->Fill(eventClass, pair.qLong);
          bf.fhSameQout[charge1][charge2]->Fill(eventClass, pair.qOut);
          bf.fhSameQside[charge1][charge2]->Fill(eventClass, pair.qSide);
        } else {
          bf.fhMixedQinv[charge1][charge2]->Fill(eventClass, pair.qInv);
          bf.fhMixedQlong[charge1][charge2]->Fill(eventClass, pair.qLong);
          bf.fhMixedQout[charge1][charge2]->Fill(eventClass, pair.qOut);
          bf.fhMixedQside[charge1][charge2]->Fill(eventClass, pair.qSide);
        }
      }
    }
  }

  /// Same-event balance-function inputs of one collision: event and kaon counters, factorial moments, pairs
  template <typename TParticles>
  void fillBfSameEvent(BfHistos& bf, int eventClass, TParticles const& particles)
  {
    if (eventClass < 0) {
      bf.fhEventsOutsideClasses->Fill(1.);
      return;
    }
    std::array<int, NCharges> nKaons{0, 0};
    for (auto const& particle : particles) {
      const int charge = particle.trackacceptedid();
      if (charge < NCharges) {
        nKaons[charge]++;
      }
    }
    bf.fhEventsPerClass->Fill(eventClass);
    for (int charge1 = 0; charge1 < NCharges; ++charge1) {
      bf.fhKaonsPerClass->Fill(eventClass, charge1, nKaons[charge1]);
      for (int charge2 = 0; charge2 < NCharges; ++charge2) {
        const int pairsInEvent = nKaons[charge1] * (nKaons[charge2] - (charge1 == charge2 ? 1 : 0));
        bf.fhFactorialPerClass->Fill(eventClass, charge1 * NCharges + charge2, pairsInEvent);
      }
    }
    fillBfPairs<true>(bf, eventClass, particles, particles);
  }

  /// Decay mother of an MC particle, looked up in the particle table given (one MC collision or the whole table).
  /// Only hadronic, non-beam mothers count: in the stored MC record, hadrons produced directly in string fragmentation
  /// point to a beam proton, which would otherwise make unrelated kaons share a "mother".
  template <typename TParticle, typename TAllParticles>
  static MotherInfo findMother(TParticle const& particle, TAllParticles const& allParticles)
  {
    MotherInfo info;
    const auto& motherIds = particle.mothersIds();
    if (motherIds.size() == 0 || motherIds[0] < 0) {
      return info;
    }
    const int64_t localIndex = motherIds[0] - static_cast<int64_t>(allParticles.offset());
    if (localIndex < 0 || localIndex >= static_cast<int64_t>(allParticles.size())) {
      return info; // mother outside the table given: treat as unknown
    }
    const auto mother = allParticles.rawIteratorAt(localIndex);
    const int genStatus = mother.getGenStatusCode();
    const bool isBeam = (genStatus >= KBeamStatusCodeMin && genStatus <= KBeamStatusCodeMax) || mother.getHepMCStatusCode() == KHepMCBeamStatusCode;
    if (isBeam || std::abs(mother.pdgCode()) <= KPdgHadronThreshold) {
      return info; // beam particle or parton/string entry: not a decay mother, the particle counts as having no mother
    }
    info.index = motherIds[0];
    info.pdgCode = mother.pdgCode();
    return info;
  }

  /// Truth of one MC particle
  template <typename TMcParticle, typename TAllParticles>
  static TruthInfo truthOf(TMcParticle const& mcParticle, TAllParticles const& allParticles)
  {
    TruthInfo truth;
    truth.hasMc = true;
    truth.mcIndex = mcParticle.globalIndex();
    truth.pdgCode = mcParticle.pdgCode();
    truth.pt = mcParticle.pt();
    truth.eta = mcParticle.eta();
    truth.phi = mcParticle.phi();
    truth.mother = findMother(mcParticle, allParticles);
    return truth;
  }

  /// Fills singles and ordered pair histograms of one level for one collision; truth[i] belongs to the i-th particle
  template <typename TCollision, typename TParticles>
  void fillCollision(PairHistos& histos, TCollision const& collision, TParticles const& particles, std::vector<TruthInfo> const& truth)
  {
    const float centMult = collision.centmult();
    const double massKaon = constants::physics::MassKaonCharged;
    histos.fhEvents->Fill(centMult, collision.posZ());

    std::size_t position1 = 0;
    for (auto const& particle1 : particles) {
      const std::size_t index1 = position1++;
      const int charge1 = particle1.trackacceptedid();
      if (charge1 >= NCharges) {
        histos.fhIdCheck->Fill(kIdCheckNotKaonId);
        continue;
      }
      const TruthInfo& truth1 = truth[index1];
      const int expectedPdg1 = (charge1 == 0) ? kKPlus : kKMinus;
      const bool trueKaon1 = truth1.hasMc && truth1.pdgCode == expectedPdg1;
      histos.fhIdCheck->Fill(!truth1.hasMc ? kIdCheckNoMcTruth : (trueKaon1 ? kIdCheckGood : kIdCheckWrongPdg));
      histos.fhPtVsCent[charge1]->Fill(centMult, particle1.pt());
      histos.fhEta[charge1]->Fill(particle1.eta());
      histos.fhPhi[charge1]->Fill(particle1.phi());
      if (truth1.mother.pdgCode == o2::constants::physics::Pdg::kPhi) {
        histos.fhPtFromPhi[charge1]->Fill(particle1.pt());
      }

      std::size_t position2 = 0;
      for (auto const& particle2 : particles) {
        const std::size_t index2 = position2++;
        const int charge2 = particle2.trackacceptedid();
        if (charge2 >= NCharges || particle1.globalIndex() == particle2.globalIndex()) {
          continue;
        }
        const TruthInfo& truth2 = truth[index2];
        const auto pair = analysis::qspace::computePairQ(particle1.pt(), particle1.eta(), particle1.phi(), massKaon,
                                                         particle2.pt(), particle2.eta(), particle2.phi(), massKaon);
        histos.fhQinv[charge1][charge2]->Fill(centMult, pair.qInv);
        histos.fhQlong[charge1][charge2]->Fill(centMult, pair.qLong);
        histos.fhQout[charge1][charge2]->Fill(centMult, pair.qOut);
        histos.fhQside[charge1][charge2]->Fill(centMult, pair.qSide);
        histos.fhMass[charge1][charge2]->Fill(pair.mass);
        const double qInvSquared = pair.qInv * pair.qInv;
        histos.fhSumRuleResidual[charge1][charge2]->Fill(pair.qLong * pair.qLong + pair.qOut * pair.qOut + pair.qSide * pair.qSide - qInvSquared);
        histos.fhMassRuleResidual[charge1][charge2]->Fill(pair.mass * pair.mass - 4. * massKaon * massKaon - qInvSquared);

        // origin: both truths known, two different MC particles (two tracks of one split particle are not a decay pair), same decay mother
        const bool twoMcParticles = truth1.hasMc && truth2.hasMc && truth1.mcIndex != truth2.mcIndex;
        const bool commonMother = twoMcParticles && truth1.mother.index >= 0 && truth1.mother.index == truth2.mother.index;
        int origin = kPairNoCommonMother;
        if (commonMother) {
          origin = (truth1.mother.pdgCode == o2::constants::physics::Pdg::kPhi) ? kPairPhi : kPairOtherCommonMother;
          histos.fhCommonMotherPdg->Fill(std::abs(truth1.mother.pdgCode));
        }
        histos.fhQinvOrigin[origin][charge1][charge2]->Fill(centMult, pair.qInv);
        histos.fhQlongOrigin[origin][charge1][charge2]->Fill(centMult, pair.qLong);
        histos.fhQoutOrigin[origin][charge1][charge2]->Fill(centMult, pair.qOut);
        histos.fhQsideOrigin[origin][charge1][charge2]->Fill(centMult, pair.qSide);

        if (histos.withResolution) {
          const int expectedPdg2 = (charge2 == 0) ? kKPlus : kKMinus;
          int status = kPairBothTrueKaons;
          if (!truth1.hasMc || !truth2.hasMc) {
            status = kPairUnmatched;
          } else if (truth1.mcIndex == truth2.mcIndex) {
            status = kPairSameMcParticle;
          } else if (!trueKaon1 || truth2.pdgCode != expectedPdg2) {
            status = kPairMisidentified;
          }
          histos.fhPairTruth[charge1][charge2]->Fill(status);
          if (status == kPairBothTrueKaons) {
            const auto pairTrue = analysis::qspace::computePairQ(truth1.pt, truth1.eta, truth1.phi, massKaon,
                                                                 truth2.pt, truth2.eta, truth2.phi, massKaon);
            histos.fhResQinv[charge1][charge2]->Fill(pairTrue.qInv, pair.qInv - pairTrue.qInv);
            histos.fhResQlong[charge1][charge2]->Fill(pairTrue.qLong, pair.qLong - pairTrue.qLong);
            histos.fhResQout[charge1][charge2]->Fill(pairTrue.qOut, pair.qOut - pairTrue.qOut);
            histos.fhResQside[charge1][charge2]->Fill(pairTrue.qSide, pair.qSide - pairTrue.qSide);
          }
        }
      }
    }
  }

  void processGenLevelNotStored(soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>::iterator const& collision,
                                soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>> const& particles,
                                soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const& allParticles)
  {
    std::vector<TruthInfo> truth;
    truth.reserve(particles.size());
    for (auto const& particle : particles) {
      truth.push_back(truthOf(particle, allParticles));
    }
    fillCollision(genHistos, collision, particles, truth);
    fillBfSameEvent(genHistos.bf, binningGen.getBin(std::make_tuple(collision.posZ(), collision.centmult())), particles);
  }
  PROCESS_SWITCH(QSpaceBalanceFunction, processGenLevelNotStored, "Generator level, filter tables made in the same workflow", true);

  void processRecLevelMcNotStored(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                  soa::Filtered<soa::Join<aod::Tracks, aod::DptDptCFTracksInfo, aod::McTrackLabels>> const& tracks,
                                  aod::McParticles const& mcParticles)
  {
    std::vector<TruthInfo> truth;
    truth.reserve(tracks.size());
    for (auto const& track : tracks) {
      if (track.has_mcParticle()) {
        truth.push_back(truthOf(track.template mcParticle_as<aod::McParticles>(), mcParticles));
      } else {
        truth.emplace_back();
      }
    }
    fillCollision(recHistos, collision, tracks, truth);
    fillBfSameEvent(recHistos.bf, binningRec.getBin(std::make_tuple(collision.posZ(), collision.centmult())), tracks);
  }
  PROCESS_SWITCH(QSpaceBalanceFunction, processRecLevelMcNotStored, "Detector level of MC (tracks with MC labels), filter tables made in the same workflow", false);

  /// Mixed events: SameKindPair pairs each accepted collision with the next cfgMixDepth collisions of the same class
  /// (unordered, without repetition; collisions without kaons are paired too, so the event-pair count is unbiased)
  template <typename TPairs, typename TBinning>
  void mixEvents(BfHistos& bf, TPairs& eventPairs, TBinning const& binning)
  {
    for (auto const& [collision1, tracks1, collision2, tracks2] : eventPairs) {
      const int eventClass = binning.getBin(std::make_tuple(collision1.posZ(), collision1.centmult()));
      if (eventClass < 0) {
        continue;
      }
      bf.fhMixedEventPairsPerClass->Fill(eventClass);
      fillBfPairs<false>(bf, eventClass, tracks1, tracks2);
      fillBfPairs<false>(bf, eventClass, tracks2, tracks1);
    }
  }

  void processGenLevelMixedNotStored(GenCollisions const& collisions, GenParticles const& particles)
  {
    auto particlesTuple = std::make_tuple(particles);
    SameKindPair<GenCollisions, GenParticles, BinningZMultGen> eventPairs{binningGen, cfgMixDepth.value, -1, collisions, particlesTuple, &cache};
    mixEvents(genHistos.bf, eventPairs, binningGen);
  }
  PROCESS_SWITCH(QSpaceBalanceFunction, processGenLevelMixedNotStored, "Generator level mixed events, filter tables made in the same workflow", false);

  void processRecLevelMixedNotStored(RecCollisions const& collisions, RecTracks const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<RecCollisions, RecTracks, BinningZMultRec> eventPairs{binningRec, cfgMixDepth.value, -1, collisions, tracksTuple, &cache};
    mixEvents(recHistos.bf, eventPairs, binningRec);
  }
  PROCESS_SWITCH(QSpaceBalanceFunction, processRecLevelMixedNotStored, "Detector level mixed events, filter tables made in the same workflow", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<QSpaceBalanceFunction>(cfgc)};
}
