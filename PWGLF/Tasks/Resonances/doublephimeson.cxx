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
///
/// \brief this is a starting point for the Resonances tutorial
/// \author sourav kundu
/// \since 02/11/2023

#include <Framework/Configurable.h>
#include <TLorentzVector.h>
#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <fairlogger/Logger.h>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/StepTHn.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/ReducedDoublePhiTables.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct doublephimeson {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<bool> fillDeltaMass{"fillDeltaMass", 1, "Fill Delta Mass"};
  Configurable<bool> fillRotation{"fillRotation", 1, "Fill rotation"};
  Configurable<int> strategyPID1{"strategyPID1", 0, "PID strategy 1"};
  Configurable<int> strategyPID2{"strategyPID2", 0, "PID strategy 2"};
  Configurable<float> minPhiMass{"minPhiMass", 1.01, "Minimum phi mass"};
  Configurable<float> maxPhiMass{"maxPhiMass", 1.03, "Maximum phi mass"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selection"};
  Configurable<float> cutNsigmaTPC{"cutNsigmaTPC", 2.5, "nsigma cut TPC"};
  Configurable<float> cutNsigmaTOF{"cutNsigmaTOF", 3.0, "nsigma cut TOF"};
  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {10, -10, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0, 20.0, 40.0, 60.0, 80.0, 500.0}, "Mixing bins - number of contributor"};

  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {1500, 2.0, 3.5}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisInvMassPhi{"configThnAxisInvMassPhi", {20, 1.01, 1.03}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisDaugherPt{"configThnAxisDaugherPt", {25, 0.0, 50.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {40, 0.0, 20.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisKstar{"configThnAxisKstar", {200, 0.0, 2.0}, "#it{k}^{*} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisDeltaR{"configThnAxisDeltaR", {200, 0.0, 2.0}, "#it{k}^{*} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCosTheta{"configThnAxisCosTheta", {5, 0.0, 1.0}, "cos #theta{*}"};
  // ConfigurableAxis configThnAxisPhiMult{"configThnAxisPhiMult", {10, 0.5, 10.5}, "#Phi Multiplicity"};

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    // register histograms
    histos.add("hnsigmaTPCKaonPlus", "hnsigmaTPCKaonPlus", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCKaonMinus", "hnsigmaTPCKaonMinus", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hPhid1Mass", "hPhid1Mass", kTH2F, {{40, 1.0, 1.04f}, {100, 0.0f, 10.0f}});
    histos.add("hPhid2Mass", "hPhid2Mass", kTH2F, {{40, 1.0, 1.04f}, {100, 0.0f, 10.0f}});

    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassPhi{configThnAxisInvMassPhi, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisDaughterPt{configThnAxisDaugherPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisKstar{configThnAxisKstar, "#it{k}^{*} (GeV/#it{c})"};
    const AxisSpec thnAxisDeltaR{configThnAxisDeltaR, "#Delta R)"};
    const AxisSpec thnAxisCosTheta{configThnAxisCosTheta, "cos #theta{*}"};
    // const AxisSpec thnAxisPhiMult{configThnAxisPhiMult, "#Phi Multiplicity)"};
    histos.add("SEMass", "SEMass", HistType::kTHnSparseF, {thnAxisInvMassPhi, thnAxisInvMassPhi, thnAxisDaughterPt, thnAxisDaughterPt});
    histos.add("SEMassUnlike", "SEMassUnlike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisKstar, thnAxisCosTheta, thnAxisDeltaR});
    histos.add("SEMassRot", "SEMassRot", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisKstar, thnAxisCosTheta, thnAxisDeltaR});
    histos.add("MEMassUnlike", "MEMassUnlike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisKstar, thnAxisCosTheta, thnAxisDeltaR});
  }

  // get kstar
  TLorentzVector trackSum, PartOneCMS, PartTwoCMS, trackRelK;
  float getkstar(const TLorentzVector part1,
                 const TLorentzVector part2)
  {
    // const TLorentzVector trackSum = part1 + part2;
    trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    // TLorentzVector PartOneCMS(part1);
    // TLorentzVector PartTwoCMS(part2);
    PartOneCMS.SetXYZM(part1.Px(), part1.Py(), part1.Pz(), part1.M());
    PartTwoCMS.SetXYZM(part2.Px(), part2.Py(), part2.Pz(), part2.M());
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    // const TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;
    trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }

  // get cosTheta
  TLorentzVector daughterCMS;
  ROOT::Math::XYZVector threeVecDauCM, threeVecMother;
  float getCosTheta(const TLorentzVector mother,
                    const TLorentzVector daughter)
  {
    threeVecMother = mother.Vect();
    const float beta = mother.Beta();
    const float betax = beta * std::cos(mother.Phi()) * std::sin(mother.Theta());
    const float betay = beta * std::sin(mother.Phi()) * std::sin(mother.Theta());
    const float betaz = beta * std::cos(mother.Theta());
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    daughterCMS = boostPRF(daughter);
    threeVecDauCM = daughterCMS.Vect();
    float cosThetaStar = TMath::Abs(threeVecDauCM.Dot(threeVecMother) / std::sqrt(threeVecMother.Mag2()) / std::sqrt(threeVecDauCM.Mag2()));
    return cosThetaStar;
  }

  bool selectionPID(float nsigmaTPC, float nsigmaTOF, int TOFHit, int PIDStrategy, float ptcand)
  {
    if (PIDStrategy == 0) {
      if (TOFHit != 1) {
        if (TMath::Abs(nsigmaTPC) < cutNsigmaTPC) {
          return true;
        }
      }
      if (TOFHit == 1) {
        if (TMath::Abs(nsigmaTOF) < cutNsigmaTOF) {
          return true;
        }
      }
    }
    if (PIDStrategy == 1) {
      if (ptcand < 0.5) {
        if (TOFHit != 1 && TMath::Abs(nsigmaTPC) < cutNsigmaTPC) {
          return true;
        }
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cutNsigmaTOF) {
          return true;
        }
      }
      if (ptcand >= 0.5) {
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cutNsigmaTOF) {
          return true;
        }
      }
    }
    if (PIDStrategy == 2) {
      if (ptcand < 0.5) {
        if (TOFHit != 1 && TMath::Abs(nsigmaTPC) < cutNsigmaTPC) {
          return true;
        }
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cutNsigmaTOF) {
          return true;
        }
      }
      if (ptcand >= 0.5 && ptcand < 1.2) {
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cutNsigmaTOF) {
          return true;
        }
        if (TOFHit != 1 && nsigmaTPC > -1.5 && nsigmaTPC < cutNsigmaTPC) {
          return true;
        }
      }
      if (ptcand >= 1.2) {
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cutNsigmaTOF) {
          return true;
        }
        if (TOFHit != 1 && TMath::Abs(nsigmaTPC) < cutNsigmaTPC) {
          return true;
        }
      }
    }
    if (PIDStrategy == 3) {
      if (ptcand < 0.5) {
        if (TOFHit != 1 && TMath::Abs(nsigmaTPC) < cutNsigmaTPC) {
          return true;
        }
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cutNsigmaTOF) {
          return true;
        }
      }
      if (ptcand >= 0.5 && ptcand < 1.2) {
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cutNsigmaTOF) {
          return true;
        }
      }
      if (ptcand >= 1.2) {
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cutNsigmaTOF) {
          return true;
        }
        if (TOFHit != 1 && TMath::Abs(nsigmaTPC) < cutNsigmaTPC) {
          return true;
        }
      }
    }
    return false;
  }

  TLorentzVector exotic, Phid1, Phid2;
  TLorentzVector exoticRot, Phid1Rot;
  void process(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    if (additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }
    for (auto phitrackd1 : phitracks) {
      if (phitrackd1.phiMass() < minPhiMass || phitrackd1.phiMass() > maxPhiMass) {
        continue;
      }

      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());

      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), strategyPID1, kaonminusd1pt)) {
        continue;
      }
      // LOGF(info, "pass TOF hit: (%d, %d)",  phitrackd1.phid1TOFHit(),  phitrackd1.phid2TOFHit());
      histos.fill(HIST("hnsigmaTPCKaonPlus"), phitrackd1.phid1TPC(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), phitrackd1.phid2TPC(), kaonminusd1pt);
      Phid1.SetXYZM(phitrackd1.phiPx(), phitrackd1.phiPy(), phitrackd1.phiPz(), phitrackd1.phiMass());
      histos.fill(HIST("hPhid1Mass"), Phid1.M(), Phid1.Pt());
      auto phid1id = phitrackd1.index();
      for (auto phitrackd2 : phitracks) {
        auto phid2id = phitrackd2.index();
        if (phid2id <= phid1id) {
          continue;
        }
        if (phitrackd2.phiMass() < minPhiMass || phitrackd2.phiMass() > maxPhiMass) {
          continue;
        }
        auto kaonplusd2pt = TMath::Sqrt(phitrackd2.phid1Px() * phitrackd2.phid1Px() + phitrackd2.phid1Py() * phitrackd2.phid1Py());
        auto kaonminusd2pt = TMath::Sqrt(phitrackd2.phid2Px() * phitrackd2.phid2Px() + phitrackd2.phid2Py() * phitrackd2.phid2Py());

        if (!selectionPID(phitrackd2.phid1TPC(), phitrackd2.phid1TOF(), phitrackd2.phid1TOFHit(), strategyPID2, kaonplusd2pt)) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid2TPC(), phitrackd2.phid2TOF(), phitrackd2.phid2TOFHit(), strategyPID2, kaonminusd2pt)) {
          continue;
        }
        if (phitrackd1.phid1Index() == phitrackd2.phid1Index()) {
          continue;
        }
        if (phitrackd1.phid2Index() == phitrackd2.phid2Index()) {
          continue;
        }
        Phid2.SetXYZM(phitrackd2.phiPx(), phitrackd2.phiPy(), phitrackd2.phiPz(), phitrackd2.phiMass());
        exotic = Phid1 + Phid2;
        auto cosThetaStar = getCosTheta(exotic, Phid1);
        auto kstar = getkstar(Phid1, Phid2);
        auto deltaR = TMath::Sqrt(TMath::Power(Phid1.Phi() - Phid2.Phi(), 2.0) + TMath::Power(Phid1.Eta() - Phid2.Eta(), 2.0));
        if (!fillDeltaMass) {
          histos.fill(HIST("SEMassUnlike"), exotic.M(), exotic.Pt(), kstar, cosThetaStar, deltaR);
        }
        if (fillDeltaMass) {
          histos.fill(HIST("SEMassUnlike"), exotic.M() - Phid1.M(), exotic.Pt(), kstar, cosThetaStar, deltaR);
        }
        histos.fill(HIST("SEMass"), Phid1.M(), Phid2.M(), Phid1.Pt(), Phid2.Pt());
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < 5; nrotbkg++) {
            auto anglestart = 5.0 * TMath::Pi() / 6.0;
            auto angleend = 7.0 * TMath::Pi() / 6.0;
            auto anglestep = (angleend - anglestart) / (1.0 * (9.0 - 1.0));
            auto rotangle = anglestart + nrotbkg * anglestep;
            auto rotd1px = Phid1.Px() * std::cos(rotangle) - Phid1.Py() * std::sin(rotangle);
            auto rotd1py = Phid1.Px() * std::sin(rotangle) + Phid1.Py() * std::cos(rotangle);
            Phid1Rot.SetXYZM(rotd1px, rotd1py, Phid1.Pz(), Phid1.M());
            exoticRot = Phid1Rot + Phid2;
            auto cosThetaStar_rot = getCosTheta(exoticRot, Phid1Rot);
            auto kstar_rot = getkstar(Phid1Rot, Phid2);
            auto deltaR_rot = TMath::Sqrt(TMath::Power(Phid1Rot.Phi() - Phid2.Phi(), 2.0) + TMath::Power(Phid1Rot.Eta() - Phid2.Eta(), 2.0));
            if (!fillDeltaMass) {
              histos.fill(HIST("SEMassRot"), exoticRot.M(), exoticRot.Pt(), kstar_rot, cosThetaStar_rot, deltaR_rot);
            }
            if (fillDeltaMass) {
              histos.fill(HIST("SEMassRot"), exoticRot.M() - Phid1Rot.M(), exoticRot.Pt(), kstar_rot, cosThetaStar_rot, deltaR_rot);
            }
          }
        }
      }
    }
  }

  SliceCache cache;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  void processMixedEvent(aod::RedPhiEvents& collisions, aod::PhiTracks& phitracks)
  {
    auto tracksTuple = std::make_tuple(phitracks);
    BinningTypeVertexContributor binningOnPositions{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::RedPhiEvents, aod::PhiTracks, BinningTypeVertexContributor> pair{binningOnPositions, nEvtMixing, -1, collisions, tracksTuple, &cache};

    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (collision1.index() == collision2.index()) {
        continue;
      }
      if (additionalEvsel && (collision1.numPos() < 2 || collision1.numNeg() < 2)) {
        continue;
      }
      if (additionalEvsel && (collision2.numPos() < 2 || collision2.numNeg() < 2)) {
        continue;
      }
      for (auto& [phitrackd1, phitrackd2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (phitrackd1.phiMass() < minPhiMass || phitrackd1.phiMass() > maxPhiMass) {
          continue;
        }
        if (phitrackd2.phiMass() < minPhiMass || phitrackd2.phiMass() > maxPhiMass) {
          continue;
        }
        auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
        auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
        auto kaonplusd2pt = TMath::Sqrt(phitrackd2.phid1Px() * phitrackd2.phid1Px() + phitrackd2.phid1Py() * phitrackd2.phid1Py());
        auto kaonminusd2pt = TMath::Sqrt(phitrackd2.phid2Px() * phitrackd2.phid2Px() + phitrackd2.phid2Py() * phitrackd2.phid2Py());
        if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), strategyPID1, kaonplusd1pt)) {
          continue;
        }
        if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), strategyPID1, kaonminusd1pt)) {
          continue;
        }
        Phid1.SetXYZM(phitrackd1.phiPx(), phitrackd1.phiPy(), phitrackd1.phiPz(), phitrackd1.phiMass());
        if (!selectionPID(phitrackd2.phid1TPC(), phitrackd2.phid1TOF(), phitrackd2.phid1TOFHit(), strategyPID2, kaonplusd2pt)) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid2TPC(), phitrackd2.phid2TOF(), phitrackd2.phid2TOFHit(), strategyPID2, kaonminusd2pt)) {
          continue;
        }
        Phid2.SetXYZM(phitrackd2.phiPx(), phitrackd2.phiPy(), phitrackd2.phiPz(), phitrackd2.phiMass());
        exotic = Phid1 + Phid2;
        auto cosThetaStar = getCosTheta(exotic, Phid1);
        auto kstar = getkstar(Phid1, Phid2);
        auto deltaR = TMath::Sqrt(TMath::Power(Phid1.Phi() - Phid2.Phi(), 2.0) + TMath::Power(Phid1.Eta() - Phid2.Eta(), 2.0));
        if (!fillDeltaMass) {
          histos.fill(HIST("MEMassUnlike"), exotic.M(), exotic.Pt(), kstar, cosThetaStar, deltaR);
        }
        if (fillDeltaMass) {
          histos.fill(HIST("MEMassUnlike"), exotic.M() - Phid1.M(), exotic.Pt(), kstar, cosThetaStar, deltaR);
        }
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processMixedEvent, "Process EventMixing for combinatorial background", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<doublephimeson>(cfgc)}; }
