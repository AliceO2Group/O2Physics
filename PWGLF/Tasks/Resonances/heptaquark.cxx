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
/// \author junlee kim

#include "PWGLF/DataModel/ReducedHeptaQuarkTables.h"

#include "Common/Core/trackUtilities.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>

#include <fairlogger/Logger.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct heptaquark {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<bool> cfgRotBkg{"cfgRotBkg", true, "flag to construct rotational backgrounds"};
  Configurable<int> cfgNRotBkg{"cfgNRotBkg", 10, "the number of rotational backgrounds"};

  Configurable<int> cfgPIDStrategy{"cfgPIDStrategy", 3, "PID strategy 1"};
  Configurable<float> cfgPIDPrPi{"cfgPIDPrPi", 3, "PID selection for proton and pion"};

  Configurable<float> cfgMinPhiMass{"cfgMinPhiMass", 1.01, "Minimum phi mass"};
  Configurable<float> cfgMaxPhiMass{"cfgMaxPhiMass", 1.03, "Maximum phi mass"};

  Configurable<float> cfgMinLambdaMass{"cfgMinLambdaMass", 1.1, "Minimum lambda mass"};
  Configurable<float> cfgMaxLambdaMass{"cfgMaxLambdaMass", 1.13, "Maximum lambda mass"};

  Configurable<float> cfgNsigmaTPC{"cfgNsigmaTPC", 2.5, "nsigma cut TPC"};
  Configurable<float> cfgNsigmaTOF{"cfgNsigmaTOF", 3.0, "nsigma cut TOF"};

  Configurable<bool> cfgSelectHQ{"cfgSelectHQ", true, "switch to select HQ"};

  Configurable<float> cfgMinPhiPt{"cfgMinPhiPt", 0.2, "Minimum phi pt"};
  Configurable<float> cfgMinLambdaPt{"cfgMinLambdaPt", 0.5, "Minimum lambda pt"};

  Configurable<float> cfgSoftFraction{"cfgSoftFraction", 0.01, "Minimum allowed softest fraction"};
  Configurable<float> cfgCollinear{"cfgCollinear", 0.98, "Maximum allowed collinear selection"};
  Configurable<float> cfgCosPoint{"cfgCosPoint", 0.95, "Minimum pointing angle selection"};

  ConfigurableAxis massAxis{"massAxis", {600, 2.8, 3.4}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "Transverse momentum bins"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0, 20, 50, 100}, "Centrality interval"};
  ConfigurableAxis massPPAxis{"massPPAxis", {100, 3.0, 8.0}, "Mass square of phi phi"};
  ConfigurableAxis massPLAxis{"massPLAxis", {100, 3.0, 8.0}, "Mass square of phi lambda"};

  void init(o2::framework::InitContext&)
  {
    histos.add("hnsigmaTPCPi", "hnsigmaTPCPi", kTH2F, {{1000, -7.0, 7.0f}, {100, 0.0f, 10.0f}});

    histos.add("hnsigmaTPCKa", "hnsigmaTPCKa", kTH2F, {{1000, -7.0, 7.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTOFKa", "hnsigmaTOFKa", kTH2F, {{1000, -7.0, 7.0f}, {100, 0.0f, 10.0f}});

    histos.add("hnsigmaTPCPr", "hnsigmaTPCPr", kTH2F, {{1000, -7.0, 7.0f}, {100, 0.0f, 10.0f}});

    histos.add("hPhid1Mass", "hPhid1Mass", kTH2F, {{40, 1.0, 1.04f}, {100, 0.0f, 10.0f}});
    histos.add("hPhid2Mass", "hPhid2Mass", kTH2F, {{40, 1.0, 1.04f}, {100, 0.0f, 10.0f}});
    histos.add("hLambdaMass", "hLambdaMass", kTH2F, {{140, 1.095, 1.135}, {100, 0.0f, 10.0f}});

    histos.add("h_InvMass_same", "h_InvMass_same", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("h_InvMass_rotPhi", "h_InvMass_rotPhi", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("h_InvMass_rotLambda", "h_InvMass_rotLambda", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("h_InvMass_rotPhiLambda", "h_InvMass_rotPhiLambda", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});

    histos.add("hDalitz", "hDalitz", {HistType::kTHnSparseF, {massPPAxis, massPLAxis, massAxis, ptAxis, {2, -0.5f, 1.5f}, centAxis}});
    histos.add("hDalitzRot", "hDalitzRot", {HistType::kTHnSparseF, {massPPAxis, massPLAxis, massAxis, ptAxis, {2, -0.5f, 1.5f}, centAxis}});
  }

  template <typename HQRow>
  static inline TLorentzVector makeP4FromHQRow(HQRow const& hq)
  {
    const double px = hq.hqPx();
    const double py = hq.hqPy();
    const double pz = hq.hqPz();
    const double m  = hq.hqMass();
    TLorentzVector v;
    v.SetXYZM(px, py, pz, m);
    return v;
  }

  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;
  double massKa = o2::constants::physics::MassKPlus;

  TRandom* rn = new TRandom();

  bool selectionPIDKaon(float nsigmaTPC, float nsigmaTOF, int TOFHit, int PIDStrategy, float ptcand)
  {
    if (PIDStrategy == 0) {
      if (TOFHit != 1) {
        if (TMath::Abs(nsigmaTPC) < cfgNsigmaTPC) {
          return true;
        }
      }
      if (TOFHit == 1) {
        if (TMath::Abs(nsigmaTOF) < cfgNsigmaTOF) {
          return true;
        }
      }
    }
    if (PIDStrategy == 1) {
      if (ptcand < 0.5) {
        if (TOFHit != 1 && TMath::Abs(nsigmaTPC) < cfgNsigmaTPC) {
          return true;
        }
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cfgNsigmaTOF) {
          return true;
        }
      }
      if (ptcand >= 0.5) {
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cfgNsigmaTOF) {
          return true;
        }
      }
    }
    if (PIDStrategy == 2) {
      if (ptcand < 0.5) {
        if (TOFHit != 1 && TMath::Abs(nsigmaTPC) < cfgNsigmaTPC) {
          return true;
        }
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cfgNsigmaTOF) {
          return true;
        }
      }
      if (ptcand >= 0.5 && ptcand < 1.2) {
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cfgNsigmaTOF) {
          return true;
        }
        if (TOFHit != 1 && nsigmaTPC > -1.5 && nsigmaTPC < cfgNsigmaTPC) {
          return true;
        }
      }
      if (ptcand >= 1.2) {
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cfgNsigmaTOF) {
          return true;
        }
        if (TOFHit != 1 && TMath::Abs(nsigmaTPC) < cfgNsigmaTPC) {
          return true;
        }
      }
    }
    if (PIDStrategy == 3) {
      if (ptcand < 0.5) {
        if (TOFHit != 1 && TMath::Abs(nsigmaTPC) < cfgNsigmaTPC) {
          return true;
        }
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cfgNsigmaTOF) {
          return true;
        }
      }
      if (ptcand >= 0.5 && ptcand < 1.2) {
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cfgNsigmaTOF) {
          return true;
        }
      }
      if (ptcand >= 1.2) {
        if (TOFHit == 1 && TMath::Abs(nsigmaTOF) < cfgNsigmaTOF) {
          return true;
        }
        if (TOFHit != 1 && TMath::Abs(nsigmaTPC) < cfgNsigmaTPC) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename HQRow1, typename HQRow2, typename HQRow3, typename ColRow>
  int selectHQ(HQRow1 const& hq1r, HQRow2 const& hq2r, HQRow3 const& hq3r, ColRow const& col)
  {
    int selection = 0;

    auto hq1 = makeP4FromHQRow(hq1r);
    auto hq2 = makeP4FromHQRow(hq2r);
    auto hq3 = makeP4FromHQRow(hq3r);

    if (hq1.Pt() < cfgMinPhiPt || hq2.Pt() < cfgMinPhiPt || hq3.Pt() < cfgMinLambdaPt) {
      selection += 1;
    }

    const double sumE = hq1.E() + hq2.E() + hq3.E();
    const double emin = std::min({hq1.E(), hq2.E(), hq3.E()});
    const double fmin = emin / std::max(1e-9, sumE);
    if (fmin < cfgSoftFraction) {
      selection += 2;
    }

    auto ex = hq1 + hq2 + hq3;
    TVector3 boost = -ex.BoostVector();

    auto hqphipair_boost = hq1 + hq2;
    auto hqlambda_boost  = hq3;
    hqphipair_boost.Boost(boost);
    hqlambda_boost.Boost(boost);

    const double denom = (hqlambda_boost.Vect().Mag() * hqphipair_boost.Vect().Mag());
    const double cosHel = (denom > 0.) ? (hqlambda_boost.Vect().Dot(hqphipair_boost.Vect()) / denom) : 1.0;
    if (std::abs(cosHel) > cfgCollinear) {
      selection += 4;
    }

    ROOT::Math::XYZVector rPV(col.posX(), col.posY(), col.posZ());
    ROOT::Math::XYZVector rSV(hq3r.hqx(), hq3r.hqy(), hq3r.hqz());

    ROOT::Math::XYZVector L = rSV - rPV;
    ROOT::Math::XYZVector exMom(ex.Px(), ex.Py(), ex.Pz());

    const double denom2 = (L.R() * exMom.R() + 1e-9);
    const double cosPoint = L.Dot(exMom) / denom2;
    if (cosPoint < cfgCosPoint) {
      selection += 8;
    }

    return selection;
  }

  ROOT::Math::PxPyPzMVector DauVec1, DauVec2;

  TLorentzVector exotic, HQ1, HQ2, HQ3;
  TLorentzVector exoticRot2, HQ2Rot;
  TLorentzVector exoticRot3, HQ3Rot;
  TLorentzVector exoticRot23;
  TLorentzVector HQ12, HQ13;
  TLorentzVector HQ12Rot, HQ13Rot;

  void processSameEvent(aod::RedHQEvents::iterator const& collision, aod::HQTracks const& hqtracks)
  {
    if (collision.numLambda() < 1 || collision.numPhi() < 2)
      return;

    for (auto hqtrackd1 : hqtracks) {
      if (hqtrackd1.hqId() != 333)
        continue;

      if (hqtrackd1.hqMass() < cfgMinPhiMass || hqtrackd1.hqMass() > cfgMaxPhiMass)
        continue;

      DauVec1 = ROOT::Math::PxPyPzMVector(hqtrackd1.hqd1Px(), hqtrackd1.hqd1Py(), hqtrackd1.hqd1Pz(), massKa);
      DauVec2 = ROOT::Math::PxPyPzMVector(hqtrackd1.hqd2Px(), hqtrackd1.hqd2Py(), hqtrackd1.hqd2Pz(), massKa);

      if (!selectionPIDKaon(hqtrackd1.hqd1TPC(), hqtrackd1.hqd1TOF(), hqtrackd1.hqd1TOFHit(), cfgPIDStrategy, DauVec1.Pt()))
        continue;

      if (!selectionPIDKaon(hqtrackd1.hqd2TPC(), hqtrackd1.hqd2TOF(), hqtrackd1.hqd2TOFHit(), cfgPIDStrategy, DauVec2.Pt()))
        continue;

      HQ1.SetXYZM(hqtrackd1.hqPx(), hqtrackd1.hqPy(), hqtrackd1.hqPz(), hqtrackd1.hqMass());

      histos.fill(HIST("hnsigmaTPCKa"), hqtrackd1.hqd1TPC(), DauVec1.Pt());
      histos.fill(HIST("hnsigmaTPCKa"), hqtrackd1.hqd2TPC(), DauVec2.Pt());
      if (hqtrackd1.hqd1TOFHit())
        histos.fill(HIST("hnsigmaTOFKa"), hqtrackd1.hqd1TOF(), DauVec1.Pt());
      if (hqtrackd1.hqd2TOFHit())
        histos.fill(HIST("hnsigmaTOFKa"), hqtrackd1.hqd2TOF(), DauVec2.Pt());

      auto hqd1id = hqtrackd1.index();
      histos.fill(HIST("hPhid1Mass"), HQ1.M(), HQ1.Pt());

      for (auto hqtrackd2 : hqtracks) {
        auto hqd2id = hqtrackd2.index();
        if (hqd2id <= hqd1id)
          continue;

        if (hqtrackd2.hqId() != 333)
          continue;

        if (hqtrackd2.hqMass() < cfgMinPhiMass || hqtrackd2.hqMass() > cfgMaxPhiMass)
          continue;

        DauVec1 = ROOT::Math::PxPyPzMVector(hqtrackd2.hqd1Px(), hqtrackd2.hqd1Py(), hqtrackd2.hqd1Pz(), massKa);
        DauVec2 = ROOT::Math::PxPyPzMVector(hqtrackd2.hqd2Px(), hqtrackd2.hqd2Py(), hqtrackd2.hqd2Pz(), massKa);

        if (!selectionPIDKaon(hqtrackd2.hqd1TPC(), hqtrackd2.hqd1TOF(), hqtrackd2.hqd1TOFHit(), cfgPIDStrategy, DauVec1.Pt()))
          continue;

        if (!selectionPIDKaon(hqtrackd2.hqd2TPC(), hqtrackd2.hqd2TOF(), hqtrackd2.hqd2TOFHit(), cfgPIDStrategy, DauVec2.Pt()))
          continue;

        if (hqtrackd1.hqd1Index() == hqtrackd2.hqd1Index())
          continue;

        if (hqtrackd1.hqd2Index() == hqtrackd2.hqd2Index())
          continue;

        HQ2.SetXYZM(hqtrackd2.hqPx(), hqtrackd2.hqPy(), hqtrackd2.hqPz(), hqtrackd2.hqMass());

        histos.fill(HIST("hnsigmaTPCKa"), hqtrackd2.hqd1TPC(), DauVec1.Pt());
        histos.fill(HIST("hnsigmaTPCKa"), hqtrackd2.hqd2TPC(), DauVec2.Pt());
        if (hqtrackd2.hqd1TOFHit())
          histos.fill(HIST("hnsigmaTOFKa"), hqtrackd2.hqd1TOF(), DauVec1.Pt());
        if (hqtrackd2.hqd2TOFHit())
          histos.fill(HIST("hnsigmaTOFKa"), hqtrackd2.hqd2TOF(), DauVec2.Pt());
        histos.fill(HIST("hPhid2Mass"), HQ2.M(), HQ2.Pt());

        for (auto hqtrackd3 : hqtracks) {
          if (std::abs(hqtrackd3.hqId()) != 3122)
            continue;

          if (hqtrackd3.hqMass() < cfgMinLambdaMass || hqtrackd3.hqMass() > cfgMaxLambdaMass)
            continue;

          int isLambda = static_cast<int>(hqtrackd3.hqId() < 0);

          if (hqtrackd3.hqId() > 0) {
            DauVec1 = ROOT::Math::PxPyPzMVector(hqtrackd3.hqd1Px(), hqtrackd3.hqd1Py(), hqtrackd3.hqd1Pz(), massPr);
            DauVec2 = ROOT::Math::PxPyPzMVector(hqtrackd3.hqd2Px(), hqtrackd3.hqd2Py(), hqtrackd3.hqd2Pz(), massPi);
          } else if (hqtrackd3.hqId() < 0) {
            DauVec1 = ROOT::Math::PxPyPzMVector(hqtrackd3.hqd1Px(), hqtrackd3.hqd1Py(), hqtrackd3.hqd1Pz(), massPi);
            DauVec2 = ROOT::Math::PxPyPzMVector(hqtrackd3.hqd2Px(), hqtrackd3.hqd2Py(), hqtrackd3.hqd2Pz(), massPr);
          }

          if (std::abs(hqtrackd3.hqd1TPC()) > cfgPIDPrPi || std::abs(hqtrackd3.hqd2TPC()) > cfgPIDPrPi)
            continue;

          if (hqtrackd1.hqd1Index() == hqtrackd3.hqd1Index())
            continue;

          if (hqtrackd1.hqd2Index() == hqtrackd3.hqd2Index())
            continue;

          if (hqtrackd2.hqd1Index() == hqtrackd3.hqd1Index())
            continue;

          if (hqtrackd2.hqd2Index() == hqtrackd3.hqd2Index())
            continue;

          HQ3.SetXYZM(hqtrackd3.hqPx(), hqtrackd3.hqPy(), hqtrackd3.hqPz(), hqtrackd3.hqMass());
          histos.fill(HIST("hLambdaMass"), HQ3.M(), HQ3.Pt());

          if (hqtrackd3.hqId() > 0) {
            histos.fill(HIST("hnsigmaTPCPr"), hqtrackd3.hqd1TPC(), DauVec1.Pt());
            histos.fill(HIST("hnsigmaTPCPi"), hqtrackd3.hqd2TPC(), DauVec2.Pt());
          } else if (hqtrackd3.hqId() < 0) {
            histos.fill(HIST("hnsigmaTPCPi"), hqtrackd3.hqd1TPC(), DauVec1.Pt());
            histos.fill(HIST("hnsigmaTPCPr"), hqtrackd3.hqd2TPC(), DauVec2.Pt());
          }

          exotic = HQ1 + HQ2 + HQ3;
          HQ12 = HQ1 + HQ2;
          HQ13 = HQ1 + HQ3;

          if (cfgSelectHQ && selectHQ(hqtrackd1, hqtrackd2, hqtrackd3, collision))
            continue;

          histos.fill(HIST("h_InvMass_same"), exotic.M(), exotic.Pt(), collision.centrality());
          histos.fill(HIST("hDalitz"), HQ12.M2(), HQ13.M2(), exotic.M(), exotic.Pt(), isLambda, collision.centrality());

          if (cfgRotBkg) {
            for (int nr = 0; nr < cfgNRotBkg; nr++) {
              auto RanPhiForD2 = rn->Uniform(o2::constants::math::PI * 5.0 / 6.0, o2::constants::math::PI * 7.0 / 6.0);
              auto RanPhiForD3 = rn->Uniform(o2::constants::math::PI * 5.0 / 6.0, o2::constants::math::PI * 7.0 / 6.0);

              RanPhiForD2 += HQ2.Phi();
              RanPhiForD3 += HQ3.Phi();

              HQ2Rot.SetXYZM(HQ2.Pt() * std::cos(RanPhiForD2), HQ2.Pt() * std::sin(RanPhiForD2), HQ2.Pz(), HQ2.M());
              HQ3Rot.SetXYZM(HQ3.Pt() * std::cos(RanPhiForD3), HQ3.Pt() * std::sin(RanPhiForD3), HQ3.Pz(), HQ3.M());

              exoticRot2 = HQ1 + HQ2Rot + HQ3;
              exoticRot3 = HQ1 + HQ2 + HQ3Rot;
              exoticRot23 = HQ1 + HQ2Rot + HQ3Rot;

              HQ12Rot = HQ1 + HQ2Rot;
              HQ13Rot = HQ1 + HQ3Rot;

              histos.fill(HIST("h_InvMass_rotPhi"), exoticRot2.M(), exoticRot2.Pt(), collision.centrality());
              histos.fill(HIST("h_InvMass_rotLambda"), exoticRot3.M(), exoticRot3.Pt(), collision.centrality());
              histos.fill(HIST("h_InvMass_rotPhiLambda"), exoticRot23.M(), exoticRot23.Pt(), collision.centrality());
              histos.fill(HIST("hDalitzRot"), HQ12Rot.M2(), HQ13Rot.M2(), exoticRot23.M(), exoticRot23.Pt(), isLambda, collision.centrality());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(heptaquark, processSameEvent, "Process same event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<heptaquark>(cfgc)}; }
