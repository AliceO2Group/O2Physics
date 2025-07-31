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

#include "PWGLF/DataModel/LFSlimHeLambda.h"

#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>

#include <memory>

namespace
{
std::shared_ptr<TH2> hInvariantMassUS[2];
std::shared_ptr<TH2> hInvariantMassLS[2];
std::shared_ptr<TH2> hRotationInvariantMassUS[2];
std::shared_ptr<TH2> hRotationInvariantMassLS[2];
std::shared_ptr<TH2> hRotationInvariantMassAntiLSeta[2];
std::shared_ptr<TH2> hInvariantMassLambda[2];
std::shared_ptr<TH2> hCosPALambda;
std::shared_ptr<TH2> hNsigmaHe3;
std::shared_ptr<TH2> hNsigmaProton;
}; // namespace

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct he3LambdaDerivedAnalysis {
  HistogramRegistry mRegistry{"He3LambdaDerivedAnalysis"};

  Configurable<int> cfgNrotations{"cfgNrotations", 7, "Number of rotations for He3 candidates"};
  Configurable<bool> cfgMirrorEta{"cfgMirrorEta", true, "Mirror eta for He3 candidates"};
  Configurable<float> cfgMinCosPA{"cfgMinCosPA", 0.99, "Minimum cosPA for Lambda candidates"};
  Configurable<float> cfgMaxNSigmaTPCHe3{"cfgMaxNSigmaTPCHe3", 2.0, "Maximum nSigmaTPC for He3 candidates"};
  Configurable<float> cfgMinLambdaPt{"cfgMinLambdaPt", 0.4, "Minimum pT for Lambda candidates"};
  Configurable<float> cfgMaxLambdaDeltaM{"cfgMaxLambdaDeltaM", 10.0e-3, "Maximum deltaM for Lambda candidates"};

  void init(InitContext const&)
  {
    constexpr double ConstituentsMass = o2::constants::physics::MassProton + o2::constants::physics::MassNeutron * 2 + o2::constants::physics::MassSigmaPlus;
    for (int i = 0; i < 2; ++i) {
      hInvariantMassUS[i] = mRegistry.add<TH2>(Form("hInvariantMassUS%i", i), "Invariant Mass", {HistType::kTH2D, {{45, 1., 10}, {100, ConstituentsMass - 0.05, ConstituentsMass + 0.05}}});
      hInvariantMassLS[i] = mRegistry.add<TH2>(Form("hInvariantMassLS%i", i), "Invariant Mass", {HistType::kTH2D, {{45, 1., 10}, {100, ConstituentsMass - 0.05, ConstituentsMass + 0.05}}});
      hRotationInvariantMassUS[i] = mRegistry.add<TH2>(Form("hRotationInvariantMassUS%i", i), "Rotation Invariant Mass", {HistType::kTH2D, {{45, 1., 10}, {100, ConstituentsMass - 0.05, ConstituentsMass + 0.05}}});
      hRotationInvariantMassLS[i] = mRegistry.add<TH2>(Form("hRotationInvariantMassLS%i", i), "Rotation Invariant Mass", {HistType::kTH2D, {{45, 1., 10}, {100, ConstituentsMass - 0.05, ConstituentsMass + 0.05}}});
      hInvariantMassLambda[i] = mRegistry.add<TH2>(Form("hInvariantMassLambda%i", i), "Invariant Mass Lambda", {HistType::kTH2D, {{50, 0., 10.}, {30, o2::constants::physics::MassLambda0 - 0.015, o2::constants::physics::MassLambda0 + 0.015}}});
      hRotationInvariantMassAntiLSeta[i] = mRegistry.add<TH2>(Form("hRotationInvariantMassAntiLSeta%i", i), "Rotation Invariant Mass Anti-Lambda", {HistType::kTH2D, {{45, 1., 10}, {100, ConstituentsMass - 0.05, ConstituentsMass + 0.05}}});
    }
    hCosPALambda = mRegistry.add<TH2>("hCosPALambda", "Cosine of Pointing Angle for Lambda", {HistType::kTH2D, {{50, 0., 10.}, {500, 0.9, 1.}}});
    hNsigmaHe3 = mRegistry.add<TH2>("hNsigmaHe3", "nSigma TPC for He3", {HistType::kTH2D, {{100, -10., 10.}, {200, -5, 5.}}});
    hNsigmaProton = mRegistry.add<TH2>("hNsigmaProton", "nSigma TPC for Proton", {HistType::kTH2D, {{100, -10., 10.}, {200, -5, 5.}}});
  }

  void processSameEvent(o2::aod::LFEvents::iterator const& collision, o2::aod::LFHe3_000 const& he3s, o2::aod::LFLambda_000 const& lambdas)
  {
    std::vector<he3Candidate> he3Candidates;
    he3Candidates.reserve(he3s.size());
    std::vector<lambdaCandidate> lambdaCandidates;
    lambdaCandidates.reserve(lambdas.size());
    for (const auto& he3 : he3s) {
      if (he3.lfEventId() != collision.globalIndex()) {
        std::cout << "He3 candidate does not match event index, skipping." << std::endl;
        return;
      }
      he3Candidate candidate;
      candidate.momentum = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(he3.pt(), he3.eta(), he3.phi(), o2::constants::physics::MassHelium3);
      candidate.nSigmaTPC = he3.nSigmaTPC();
      candidate.dcaXY = he3.dcaXY();
      candidate.dcaZ = he3.dcaZ();
      candidate.tpcNClsFound = he3.tpcNCls();
      candidate.itsNCls = he3.itsClusterSizes();
      candidate.itsClusterSizes = he3.itsClusterSizes();
      candidate.sign = he3.sign();
      hNsigmaHe3->Fill(he3.pt() * he3.sign(), he3.nSigmaTPC());
      if (std::abs(he3.nSigmaTPC()) > cfgMaxNSigmaTPCHe3) {
        continue; // Skip candidates with nSigmaTPC outside range
      }
      he3Candidates.push_back(candidate);
    }
    for (const auto& lambda : lambdas) {
      if (lambda.lfEventId() != collision.globalIndex()) {
        std::cout << "Lambda candidate does not match event index, skipping." << std::endl;
        return;
      }
      lambdaCandidate candidate;
      candidate.momentum.SetCoordinates(lambda.pt(), lambda.eta(), lambda.phi(), o2::constants::physics::MassLambda0);
      candidate.mass = lambda.mass();
      candidate.cosPA = lambda.cosPA();
      candidate.dcaV0Daughters = lambda.dcaDaughters();
      hCosPALambda->Fill(lambda.pt(), candidate.cosPA);
      // hNsigmaProton->Fill(lambda.pt() * lambda.sign(), lambda.protonNSigmaTPC());
      hInvariantMassLambda[0]->Fill(lambda.pt(), lambda.mass());
      if (candidate.cosPA < cfgMinCosPA || lambda.pt() < cfgMinLambdaPt ||
          std::abs(lambda.mass() - o2::constants::physics::MassLambda0) > cfgMaxLambdaDeltaM) {
        continue; // Skip candidates with low cosPA
      }
      hInvariantMassLambda[1]->Fill(lambda.pt(), lambda.mass());
      lambdaCandidates.push_back(candidate);
    }

    for (const auto& he3 : he3Candidates) {
      for (const auto& lambda : lambdaCandidates) {
        auto pairMomentum = lambda.momentum + he3.momentum; // Calculate invariant mass
        (he3.sign * lambda.sign > 0 ? hInvariantMassLS : hInvariantMassUS)[he3.sign > 0]->Fill(pairMomentum.Pt(), pairMomentum.M());
      }
      for (int iEta{0}; iEta <= cfgMirrorEta; ++iEta) {
        for (int iR{0}; iR <= cfgNrotations; ++iR) {
          auto he3Momentum = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(he3.momentum.Pt(), (1. - iEta * 2.) * he3.momentum.Eta(), he3.momentum.Phi() + TMath::Pi() * (0.75 + 0.5 * iR / cfgNrotations), he3.momentum.M());
          for (const auto& lambda : lambdaCandidates) {
            auto pairMomentum = lambda.momentum + he3Momentum; // Calculate invariant mass
            (he3.sign * lambda.sign > 0 ? hRotationInvariantMassLS : hRotationInvariantMassUS)[he3.sign > 0]->Fill(pairMomentum.Pt(), pairMomentum.M());
            if (he3.sign < 0 && lambda.sign < 0) {
              hRotationInvariantMassAntiLSeta[iEta]->Fill(pairMomentum.Pt(), pairMomentum.M());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(he3LambdaDerivedAnalysis, processSameEvent, "Process same event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<he3LambdaDerivedAnalysis>(cfgc)};
}
