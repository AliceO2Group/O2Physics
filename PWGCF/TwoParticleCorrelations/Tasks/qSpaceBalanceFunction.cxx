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
/// At generator level every pair is also classified by the mother particles of the two kaons:
/// phi pair (same mother, mother is a phi(1020)), other common mother, or no common mother.

#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TwoParticleCorrelations/Core/QSpaceKinematics.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
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
  kIdCheckWrongPdg
};
/// Pair origin classes from the MC mother of the two kaons
enum PairOrigin {
  kPairPhi = 0,
  kPairOtherCommonMother,
  kPairNoCommonMother,
  NPairOrigins
};
const std::array<std::string, NPairOrigins> pairOriginTags{"phi", "otherMother", "noCommonMother"};
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
} // namespace

struct QSpaceBalanceFunction {
  ConfigurableAxis cfgBinsQinv{"cfgBinsQinv", {2500, 0., 5.}, "Q_inv binning (GeV/c)"};
  ConfigurableAxis cfgBinsQComponent{"cfgBinsQComponent", {1000, -5., 5.}, "Q_long, Q_out, Q_side binning (GeV/c)"};
  ConfigurableAxis cfgBinsMass{"cfgBinsMass", {6000, 0., 6.}, "Pair invariant mass binning (GeV/c^2)"};
  ConfigurableAxis cfgBinsPt{"cfgBinsPt", {18, 0.2, 2.}, "pT binning (GeV/c)"};
  ConfigurableAxis cfgBinsEta{"cfgBinsEta", {16, -0.8, 0.8}, "Pseudorapidity binning"};
  ConfigurableAxis cfgBinsResidual{"cfgBinsResidual", {200, -1.e-10, 1.e-10}, "Residual binning (GeV^2/c^2)"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  using HistoPerCharge = std::array<std::shared_ptr<TH1>, NCharges>;
  using Histo2DPerCharge = std::array<std::shared_ptr<TH2>, NCharges>;
  using Histo2DPerPair = std::array<Histo2DPerCharge, NCharges>;
  using HistoPerPair = std::array<HistoPerCharge, NCharges>;

  std::shared_ptr<TH2> fhEvents;
  std::shared_ptr<TH1> fhIdCheck;
  Histo2DPerCharge fhPtVsCent;
  HistoPerCharge fhEta, fhPhi;
  Histo2DPerPair fhQinv, fhQlong, fhQout, fhQside;
  HistoPerPair fhMass, fhSumRuleResidual, fhMassRuleResidual;
  /// pair histograms split by pair origin (generator level only): [origin][charge1][charge2]
  std::array<Histo2DPerPair, NPairOrigins> fhQinvOrigin, fhQlongOrigin, fhQoutOrigin, fhQsideOrigin;
  HistoPerCharge fhPtFromPhi;
  std::shared_ptr<TH1> fhCommonMotherPdg;

  Filter onlyAcceptedCollisions = (aod::dptdptfilter::collisionaccepted == uint8_t(true));
  Filter onlyAcceptedTracks = (aod::dptdptfilter::trackacceptedid >= int8_t(0));

  void init(InitContext const&)
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

    fhEvents = registry.add<TH2>("events/hZVsCent", "analysed events", kTH2D, {axisCent, axisZ});
    fhIdCheck = registry.add<TH1>("singles/hIdCheck", "particle id check", kTH1D, {{3, 0.5, 3.5}});
    fhIdCheck->GetXaxis()->SetBinLabel(kIdCheckGood, "kaon, id agrees with PDG");
    fhIdCheck->GetXaxis()->SetBinLabel(kIdCheckNotKaonId, "id is not 0 or 1");
    fhIdCheck->GetXaxis()->SetBinLabel(kIdCheckWrongPdg, "id disagrees with PDG");

    for (int charge1 = 0; charge1 < NCharges; ++charge1) {
      const std::string& tag1 = chargeTags[charge1];
      fhPtVsCent[charge1] = registry.add<TH2>(("singles/hPtVsCent_" + tag1).c_str(), "kaon #it{p}_{T}", kTH2D, {axisCent, axisPt});
      fhEta[charge1] = registry.add<TH1>(("singles/hEta_" + tag1).c_str(), "kaon #eta", kTH1D, {axisEta});
      fhPhi[charge1] = registry.add<TH1>(("singles/hPhi_" + tag1).c_str(), "kaon #varphi", kTH1D, {axisPhi});
      for (int charge2 = 0; charge2 < NCharges; ++charge2) {
        const std::string tag = tag1 + chargeTags[charge2];
        fhQinv[charge1][charge2] = registry.add<TH2>(("pairs/hQinv_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQinv});
        fhQlong[charge1][charge2] = registry.add<TH2>(("pairs/hQlong_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQlong});
        fhQout[charge1][charge2] = registry.add<TH2>(("pairs/hQout_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQout});
        fhQside[charge1][charge2] = registry.add<TH2>(("pairs/hQside_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQside});
        fhMass[charge1][charge2] = registry.add<TH1>(("pairs/hMass_" + tag).c_str(), "ordered pairs", kTH1F, {axisMass});
        fhSumRuleResidual[charge1][charge2] = registry.add<TH1>(("checks/hSumRuleResidual_" + tag).c_str(), "Q_{long}^{2}+Q_{out}^{2}+Q_{side}^{2}-Q_{inv}^{2}", kTH1D, {axisResidual});
        fhMassRuleResidual[charge1][charge2] = registry.add<TH1>(("checks/hMassRuleResidual_" + tag).c_str(), "M^{2}-4m^{2}-Q_{inv}^{2}", kTH1D, {axisResidual});
      }
    }
    if (doprocessGenLevelNotStored.value) {
      for (int origin = 0; origin < NPairOrigins; ++origin) {
        for (int charge1 = 0; charge1 < NCharges; ++charge1) {
          for (int charge2 = 0; charge2 < NCharges; ++charge2) {
            const std::string tag = chargeTags[charge1] + chargeTags[charge2] + "_" + pairOriginTags[origin];
            fhQinvOrigin[origin][charge1][charge2] = registry.add<TH2>(("pairsByOrigin/hQinv_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQinv});
            fhQlongOrigin[origin][charge1][charge2] = registry.add<TH2>(("pairsByOrigin/hQlong_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQlong});
            fhQoutOrigin[origin][charge1][charge2] = registry.add<TH2>(("pairsByOrigin/hQout_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQout});
            fhQsideOrigin[origin][charge1][charge2] = registry.add<TH2>(("pairsByOrigin/hQside_" + tag).c_str(), "ordered pairs", kTH2F, {axisCent, axisQside});
          }
        }
      }
      for (int charge = 0; charge < NCharges; ++charge) {
        fhPtFromPhi[charge] = registry.add<TH1>(("singles/hPtFromPhi_" + chargeTags[charge]).c_str(), "kaons whose mother is a #phi(1020)", kTH1D, {axisPt});
      }
      fhCommonMotherPdg = registry.add<TH1>("pairsByOrigin/hCommonMotherPdg", "|PDG code| of the common mother, ordered pairs (overflow = |PDG| >= 10000)", kTH1D, {{10000, 0., 10000.}});
    }
  }

  /// Decay mother of a generated particle, looked up in the unfiltered particle table of the same MC collision.
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
      return info; // mother outside this MC collision: treat as unknown
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

  /// Fills singles and ordered pair histograms for one collision
  /// mothers is empty when the mother information is not available (detector level); otherwise mothers[i] belongs to
  /// the i-th particle of the iteration
  template <typename TCollision, typename TParticles>
  void fillCollision(TCollision const& collision, TParticles const& particles, std::vector<MotherInfo> const& mothers)
  {
    const float centMult = collision.centmult();
    const double massKaon = constants::physics::MassKaonCharged;
    const bool withOrigin = !mothers.empty();
    fhEvents->Fill(centMult, collision.posZ());

    std::size_t position1 = 0;
    for (auto const& particle1 : particles) {
      const std::size_t index1 = position1++;
      const int charge1 = particle1.trackacceptedid();
      if (charge1 >= NCharges) {
        fhIdCheck->Fill(kIdCheckNotKaonId);
        continue;
      }
      fhIdCheck->Fill((particle1.pdgCode() == (charge1 == 0 ? kKPlus : kKMinus)) ? kIdCheckGood : kIdCheckWrongPdg);
      fhPtVsCent[charge1]->Fill(centMult, particle1.pt());
      fhEta[charge1]->Fill(particle1.eta());
      fhPhi[charge1]->Fill(particle1.phi());
      if (withOrigin && mothers[index1].pdgCode == o2::constants::physics::Pdg::kPhi) {
        fhPtFromPhi[charge1]->Fill(particle1.pt());
      }

      std::size_t position2 = 0;
      for (auto const& particle2 : particles) {
        const std::size_t index2 = position2++;
        const int charge2 = particle2.trackacceptedid();
        if (charge2 >= NCharges || particle1.globalIndex() == particle2.globalIndex()) {
          continue;
        }
        const auto pair = analysis::qspace::computePairQ(particle1.pt(), particle1.eta(), particle1.phi(), massKaon,
                                                         particle2.pt(), particle2.eta(), particle2.phi(), massKaon);
        fhQinv[charge1][charge2]->Fill(centMult, pair.qInv);
        fhQlong[charge1][charge2]->Fill(centMult, pair.qLong);
        fhQout[charge1][charge2]->Fill(centMult, pair.qOut);
        fhQside[charge1][charge2]->Fill(centMult, pair.qSide);
        fhMass[charge1][charge2]->Fill(pair.mass);
        const double qInvSquared = pair.qInv * pair.qInv;
        fhSumRuleResidual[charge1][charge2]->Fill(pair.qLong * pair.qLong + pair.qOut * pair.qOut + pair.qSide * pair.qSide - qInvSquared);
        fhMassRuleResidual[charge1][charge2]->Fill(pair.mass * pair.mass - 4. * massKaon * massKaon - qInvSquared);
        if (withOrigin) {
          const MotherInfo& mother1 = mothers[index1];
          const bool commonMother = (mother1.index >= 0) && (mother1.index == mothers[index2].index);
          int origin = kPairNoCommonMother;
          if (commonMother) {
            origin = (mother1.pdgCode == o2::constants::physics::Pdg::kPhi) ? kPairPhi : kPairOtherCommonMother;
            fhCommonMotherPdg->Fill(std::abs(mother1.pdgCode));
          }
          fhQinvOrigin[origin][charge1][charge2]->Fill(centMult, pair.qInv);
          fhQlongOrigin[origin][charge1][charge2]->Fill(centMult, pair.qLong);
          fhQoutOrigin[origin][charge1][charge2]->Fill(centMult, pair.qOut);
          fhQsideOrigin[origin][charge1][charge2]->Fill(centMult, pair.qSide);
        }
      }
    }
  }

  void processGenLevelNotStored(soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>::iterator const& collision,
                                soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>> const& particles,
                                soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const& allParticles)
  {
    std::vector<MotherInfo> mothers;
    mothers.reserve(particles.size());
    for (auto const& particle : particles) {
      mothers.push_back(findMother(particle, allParticles));
    }
    fillCollision(collision, particles, mothers);
  }
  PROCESS_SWITCH(QSpaceBalanceFunction, processGenLevelNotStored, "Generator level, filter tables made in the same workflow", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<QSpaceBalanceFunction>(cfgc)};
}
