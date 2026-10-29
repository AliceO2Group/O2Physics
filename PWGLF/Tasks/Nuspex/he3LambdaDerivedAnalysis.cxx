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

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <TH2.h>
#include <TMath.h>
#include <TString.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <ostream>
#include <vector>

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
std::shared_ptr<TH2> hInvariantMassHe3K0s[2]; // 0: anti-He3, 1: He3
std::shared_ptr<TH2> hRotationInvariantMassHe3K0s[2];
std::shared_ptr<TH2> hInvariantMassK0s[2]; // Before and after selection
std::shared_ptr<TH2> hCosPAK0s;
std::shared_ptr<TH2> hCtK0s;
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

  Configurable<float> cfgMinK0sPt{"cfgMinK0sPt", 0.5f, "Minimum K0s pT"};
  Configurable<float> cfgMaxK0sPt{"cfgMaxK0sPt", 10.f, "Maximum K0s pT"};
  Configurable<float> cfgMinK0sCt{"cfgMinK0sCt", 0.f, "Minimum K0s proper decay length (cm)"};
  Configurable<float> cfgMaxK0sCt{"cfgMaxK0sCt", 20.f, "Maximum K0s proper decay length (cm)"};
  Configurable<float> cfgMinK0sCosPA{"cfgMinK0sCosPA", 0.99f, "Minimum K0s cosPA"};
  Configurable<float> cfgMaxK0sDeltaM{"cfgMaxK0sDeltaM", 10.e-3f, "Maximum K0s mass deviation (GeV/c^2)"};
  Configurable<float> cfgMaxNSigmaTPCPionK0s{"cfgMaxNSigmaTPCPionK0s", 4.f, "Maximum daughter pion TPC nSigma for K0s"};
  ConfigurableAxis cfgHe3K0sMassAxis{"cfgHe3K0sMassAxis", {200, MassHelium3 + MassK0Short, MassHelium3 + MassK0Short + 1.}, "He3-K0s pair invariant mass (GeV/c^2)"};

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
    if (doprocessSameEvent && doprocessSameEventLegacy) {
      LOG(fatal) << "Enable only one Lambda table version: processSameEvent or processSameEventLegacy";
    }
    if (doprocessK0s) {
      const AxisSpec pairMassAxis{cfgHe3K0sMassAxis, "m(He3 K0s) (GeV/c^{2})"};
      for (int i = 0; i < 2; ++i) {
        hInvariantMassHe3K0s[i] = mRegistry.add<TH2>(Form("hInvariantMassHe3K0s%i", i), "Same-event He3-K0s", {HistType::kTH2D, {{45, 1., 10.}, pairMassAxis}});
        hRotationInvariantMassHe3K0s[i] = mRegistry.add<TH2>(Form("hRotationInvariantMassHe3K0s%i", i), "Rotated He3-K0s", {HistType::kTH2D, {{45, 1., 10.}, pairMassAxis}});
        hInvariantMassK0s[i] = mRegistry.add<TH2>(Form("hInvariantMassK0s%i", i), "K0s;pT (GeV/c);mass (GeV/c^{2})", {HistType::kTH2D, {{50, 0., 10.}, {100, MassK0Short - 0.015, MassK0Short + 0.015}}});
      }
      hCosPAK0s = mRegistry.add<TH2>("hCosPAK0s", "K0s;pT (GeV/c);cosPA", {HistType::kTH2D, {{50, 0., 10.}, {500, 0.9, 1.}}});
      hCtK0s = mRegistry.add<TH2>("hCtK0s", "K0s;pT (GeV/c);ct (cm)", {HistType::kTH2D, {{50, 0., 10.}, {100, 0., 30.}}});
    }
    hCosPALambda = mRegistry.add<TH2>("hCosPALambda", "Cosine of Pointing Angle for Lambda", {HistType::kTH2D, {{50, 0., 10.}, {500, 0.9, 1.}}});
    hNsigmaHe3 = mRegistry.add<TH2>("hNsigmaHe3", "nSigma TPC for He3", {HistType::kTH2D, {{100, -10., 10.}, {200, -5, 5.}}});
    hNsigmaProton = mRegistry.add<TH2>("hNsigmaProton", "nSigma TPC for Proton", {HistType::kTH2D, {{100, -10., 10.}, {200, -5, 5.}}});
  }

  template <typename He3Table, typename LambdaTable>
  void analyseLambda(o2::aod::LFEvents::iterator const& collision, He3Table const& he3s, LambdaTable const& lambdas)
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
      candidate.sign = lambda.sign();
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
        for (int iR{0}; cfgNrotations > 0 && iR <= cfgNrotations; ++iR) {
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
  void processSameEvent(o2::aod::LFEvents::iterator const& collision, o2::aod::LFHe3_001 const& he3s, o2::aod::LFLambda_001 const& lambdas)
  {
    analyseLambda(collision, he3s, lambdas);
  }
  PROCESS_SWITCH(he3LambdaDerivedAnalysis, processSameEvent, "Process same-event He3-Lambda pairs (version 001)", true);

  void processSameEventLegacy(o2::aod::LFEvents::iterator const& collision, o2::aod::LFHe3_000 const& he3s, o2::aod::LFLambda_000 const& lambdas)
  {
    analyseLambda(collision, he3s, lambdas);
  }
  PROCESS_SWITCH(he3LambdaDerivedAnalysis, processSameEventLegacy, "Process same-event He3-Lambda pairs (version 000)", false);

  void processK0s(o2::aod::LFEvents::iterator const&, o2::aod::LFHe3_001 const& he3s, o2::aod::LFK0s const& k0ss)
  {
    std::vector<k0sCandidate> k0sCandidates;
    k0sCandidates.reserve(k0ss.size());
    for (const auto& k0s : k0ss) {
      hInvariantMassK0s[0]->Fill(k0s.pt(), k0s.mass());
      hCosPAK0s->Fill(k0s.pt(), k0s.cosPA());
      hCtK0s->Fill(k0s.pt(), k0s.ct());
      if (k0s.pt() < cfgMinK0sPt || k0s.pt() > cfgMaxK0sPt || k0s.ct() < cfgMinK0sCt || k0s.ct() > cfgMaxK0sCt ||
          k0s.cosPA() < cfgMinK0sCosPA || std::abs(k0s.mass() - MassK0Short) > cfgMaxK0sDeltaM ||
          std::abs(k0s.nSigmaTPCPosPion()) > cfgMaxNSigmaTPCPionK0s || std::abs(k0s.nSigmaTPCNegPion()) > cfgMaxNSigmaTPCPionK0s) {
        continue;
      }
      k0sCandidate candidate;
      // Use the nominal K0s mass for pair kinematics, as for Lambda.
      candidate.momentum.SetCoordinates(k0s.pt(), k0s.eta(), k0s.phi(), MassK0Short);
      k0sCandidates.push_back(candidate);
      hInvariantMassK0s[1]->Fill(k0s.pt(), k0s.mass());
    }
    for (const auto& he3 : he3s) {
      if (!doprocessSameEvent && !doprocessSameEventLegacy) {
        hNsigmaHe3->Fill(he3.pt() * he3.sign(), he3.nSigmaTPC());
      }
      if (std::abs(he3.nSigmaTPC()) > cfgMaxNSigmaTPCHe3) {
        continue;
      }
      const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> he3Momentum(he3.pt(), he3.eta(), he3.phi(), MassHelium3);
      const int signIndex = he3.sign() > 0;
      for (const auto& k0s : k0sCandidates) {
        const auto pairMomentum = he3Momentum + k0s.momentum;
        hInvariantMassHe3K0s[signIndex]->Fill(pairMomentum.Pt(), pairMomentum.M());
      }
      // Match the Lambda background's rotation angles and optional eta mirroring.
      for (int iEta = 0; iEta <= cfgMirrorEta; ++iEta) {
        for (int iR = 0; cfgNrotations > 0 && iR <= cfgNrotations; ++iR) {
          const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> rotatedHe3(he3.pt(), (1. - iEta * 2.) * he3.eta(), he3.phi() + TMath::Pi() * (0.75 + 0.5 * iR / cfgNrotations), MassHelium3);
          for (const auto& k0s : k0sCandidates) {
            const auto pairMomentum = rotatedHe3 + k0s.momentum;
            hRotationInvariantMassHe3K0s[signIndex]->Fill(pairMomentum.Pt(), pairMomentum.M());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(he3LambdaDerivedAnalysis, processK0s, "Process same-event He3-K0s pairs (version 001)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<he3LambdaDerivedAnalysis>(cfgc)};
}
