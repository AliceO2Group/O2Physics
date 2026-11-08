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
/// \author Lakshitha Fernando
///
/// Reads the tables of dpt-dpt-filter, which must run in the same workflow with only kaons enabled,
/// so that trackacceptedid is 0 for K+ and 1 for K-.
/// Pairs are ORDERED, as in dpt-dpt-correlations: every particle takes the role of particle 1 in turn.
/// Each of ++, +-, -+, -- is filled; an unordered like-sign pair therefore enters twice, an unordered
/// opposite-sign pair once in +- and once in -+. OS = (+-) + (-+) and LS = (++) + (--) are comparable.

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
#include <cstdint>
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
  }

  /// Fills singles and ordered pair histograms for one collision
  template <typename TCollision, typename TParticles>
  void fillCollision(TCollision const& collision, TParticles const& particles)
  {
    const float centMult = collision.centmult();
    const double massKaon = constants::physics::MassKaonCharged;
    fhEvents->Fill(centMult, collision.posZ());

    for (auto const& particle1 : particles) {
      const int charge1 = particle1.trackacceptedid();
      if (charge1 >= NCharges) {
        fhIdCheck->Fill(kIdCheckNotKaonId);
        continue;
      }
      fhIdCheck->Fill((particle1.pdgCode() == (charge1 == 0 ? kKPlus : kKMinus)) ? kIdCheckGood : kIdCheckWrongPdg);
      fhPtVsCent[charge1]->Fill(centMult, particle1.pt());
      fhEta[charge1]->Fill(particle1.eta());
      fhPhi[charge1]->Fill(particle1.phi());

      for (auto const& particle2 : particles) {
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
      }
    }
  }

  void processGenLevelNotStored(soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>::iterator const& collision,
                                soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>> const& particles)
  {
    fillCollision(collision, particles);
  }
  PROCESS_SWITCH(QSpaceBalanceFunction, processGenLevelNotStored, "Generator level, filter tables made in the same workflow", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<QSpaceBalanceFunction>(cfgc)};
}
