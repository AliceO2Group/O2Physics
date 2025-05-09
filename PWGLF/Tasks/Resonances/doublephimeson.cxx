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
  Configurable<int> strategyPID1{"strategyPID1", 0, "PID strategy 1"};
  Configurable<int> strategyPID2{"strategyPID2", 0, "PID strategy 2"};
  Configurable<float> minPhiMass{"minPhiMass", 1.01, "Minimum phi mass"};
  Configurable<float> maxPhiMass{"maxPhiMass", 1.03, "Maximum phi mass"};
  Configurable<float> minExoticMass{"minExoticMass", 2.4, "Minimum Exotic mass"};
  Configurable<float> maxExoticMass{"maxExoticMass", 3.2, "Maximum Exotic mass"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selection"};
  Configurable<bool> isDeep{"isDeep", true, "Store deep angle"};
  Configurable<float> cutMinNsigmaTPC{"cutMinNsigmaTPC", -2.5, "nsigma cut TPC"};
  Configurable<float> cutNsigmaTPC{"cutNsigmaTPC", 3.0, "nsigma cut TPC"};
  Configurable<float> cutNsigmaTOF{"cutNsigmaTOF", 3.0, "nsigma cut TOF"};
  Configurable<float> momTOFCut{"momTOFCut", 1.8, "minimum pT cut for madnatory TOF"};
  Configurable<float> maxKaonPt{"maxKaonPt", 100.0, "maximum kaon pt cut"};
  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 1, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {10, -10, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0, 20.0, 40.0, 60.0, 80.0, 500.0}, "Mixing bins - number of contributor"};

  // THnsparse bining
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {1500, 2.0, 3.5}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisInvMassPhi{"configThnAxisInvMassPhi", {20, 1.01, 1.03}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisInvMassDeltaPhi{"configThnAxisInvMassDeltaPhi", {80, 0.0, 0.08}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisDaugherPt{"configThnAxisDaugherPt", {25, 0.0, 50.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {40, 0.0, 20.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisKstar{"configThnAxisKstar", {200, 0.0, 2.0}, "#it{k}^{*} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisDeltaR{"configThnAxisDeltaR", {200, 0.0, 2.0}, "#it{k}^{*} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCosTheta{"configThnAxisCosTheta", {160, 0.0, 3.2}, "cos #theta{*}"};
  ConfigurableAxis configThnAxisNumPhi{"configThnAxisNumPhi", {101, -0.5, 100.5}, "cos #theta{*}"};

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    // register histograms
    histos.add("hnsigmaTPCKaonPlus", "hnsigmaTPCKaonPlus", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCKaonMinus", "hnsigmaTPCKaonMinus", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCTOFKaon", "hnsigmaTPCTOFKaon", kTH3F, {{500, -3.0, 3.0f}, {500, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hPhiMass", "hPhiMass", kTH2F, {{40, 1.0, 1.04f}, {100, 0.0f, 10.0f}});
    histos.add("hPhiMass2", "hPhiMass2", kTH2F, {{40, 1.0, 1.04f}, {40, 1.0f, 1.04f}});

    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassPhi{configThnAxisInvMassPhi, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassDeltaPhi{configThnAxisInvMassDeltaPhi, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisDeltaR{configThnAxisDeltaR, "#Delta R)"};
    const AxisSpec thnAxisCosTheta{configThnAxisCosTheta, "cos #theta"};
    const AxisSpec thnAxisNumPhi{configThnAxisNumPhi, "Number of phi meson"};

    histos.add("SEMassUnlike", "SEMassUnlike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisDeltaR, thnAxisInvMassPhi, thnAxisInvMassPhi, thnAxisNumPhi});
    // histos.add("SEMassLike", "SEMassLike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisDeltaR, thnAxisInvMassPhi, thnAxisInvMassPhi, thnAxisNumPhi});
    histos.add("MEMassUnlike", "MEMassUnlike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisDeltaR, thnAxisInvMassPhi, thnAxisInvMassPhi});
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

  float deepangle2(const ROOT::Math::PtEtaPhiMVector candidate1,
                   const ROOT::Math::PtEtaPhiMVector candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.Pt();
    pt2 = candidate2.Pt();
    pz1 = candidate1.Pz();
    pz2 = candidate2.Pz();
    p1 = candidate1.P();
    p2 = candidate2.P();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    return angle;
  }

  float deepangle(const TLorentzVector candidate1,
                  const TLorentzVector candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.Pt();
    pt2 = candidate2.Pt();
    pz1 = candidate1.Pz();
    pz2 = candidate2.Pz();
    p1 = candidate1.P();
    p2 = candidate2.P();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    return angle;
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
    // optimized TPC TOF
    if (PIDStrategy == 0) {
      if (ptcand < 0.4) {
        if (nsigmaTPC > -3.0 && nsigmaTPC < 3.0) {
          return true;
        }
      } else if (ptcand >= 0.4 && ptcand < 0.5) {
        if (nsigmaTPC > -2.0 && nsigmaTPC < 3.0) {
          return true;
        }
      } else if (ptcand >= 0.5 && ptcand < 5.0 && TOFHit == 1) {
        if (ptcand < 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        }
        if (ptcand >= 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) {
          return true;
        }
      } else if (ptcand >= 0.5 && ptcand < 5.0 && TOFHit != 1) {
        if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > -0.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 1.0 && ptcand < 1.8 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 1.8 && ptcand < 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.5) {
          return true;
        }
        if (ptcand >= 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.0) {
          return true;
        }
      } else if (ptcand >= 5.0 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
        return true;
      }
    }
    // optimized TPC TOF combined
    if (PIDStrategy == 1) {
      if (ptcand < 0.4) {
        if (nsigmaTPC > cutMinNsigmaTPC && nsigmaTPC < cutNsigmaTPC) {
          return true;
        }
      } else if (ptcand >= 0.4 && ptcand < 0.5) {
        if (nsigmaTPC > -2.0 && nsigmaTPC < cutNsigmaTPC) {
          return true;
        }
      } else if (ptcand >= 0.5 && ptcand < 5.0 && TOFHit == 1) {
        if (ptcand < 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        }
        if (ptcand >= 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) {
          return true;
        }
      } else if (ptcand >= 5.0 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
        return true;
      }
    }

    if (PIDStrategy == 2) {
      if (ptcand < 0.5) {
        if (nsigmaTPC > cutMinNsigmaTPC && nsigmaTPC < cutNsigmaTPC) {
          return true;
        }
      }
      if (ptcand >= 0.5) {
        if (TOFHit != 1 && ptcand < momTOFCut) {
          if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > -0.0 && nsigmaTPC < cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 1.0 && ptcand < 1.8 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 1.8 && ptcand < 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.5) {
            return true;
          }
          if (ptcand >= 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.0) {
            return true;
          }
        }
        if (TOFHit == 1) {
          if (TMath::Sqrt((nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) / 2.0) < cutNsigmaTOF) {
            return true;
          }
        }
      }
    }
    if (PIDStrategy == 3) {
      if (ptcand < 0.5) {
        if (nsigmaTPC > cutMinNsigmaTPC && nsigmaTPC < cutNsigmaTPC) {
          return true;
        }
      }
      if (ptcand >= 0.5) {
        if (TOFHit != 1) {
          if (nsigmaTPC > cutMinNsigmaTPC && nsigmaTPC < cutNsigmaTPC) {
            return true;
          }
        }
        if (TOFHit == 1) {
          if (TMath::Sqrt((nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) / 2.0) < cutNsigmaTOF) {
            return true;
          }
        }
      }
    }
    return false;
  }

  TLorentzVector exotic, Phid1, Phid2;
  // TLorentzVector exoticlike, Phi1kaonplus, Phi1kaonminus, Phi2kaonplus, Phi2kaonminus, Phid1like, Phid2like;
  // TLorentzVector exoticRot, Phid1Rot;

  void processSE(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    if (additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }
    int phimult = 0;
    for (auto phitrackd1 : phitracks) {
      if (phitrackd1.phiMass() < minPhiMass || phitrackd1.phiMass() > maxPhiMass) {
        continue;
      }
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
      if (kaonplusd1pt > maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), strategyPID1, kaonminusd1pt)) {
        continue;
      }
      phimult = phimult + 1;
    }
    for (auto phitrackd1 : phitracks) {
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
      if (kaonplusd1pt > maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), strategyPID1, kaonplusd1pt)) {
        continue;
      }
      histos.fill(HIST("hnsigmaTPCTOFKaon"), phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), phitrackd1.phid1TPC(), kaonplusd1pt);
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), strategyPID1, kaonminusd1pt)) {
        continue;
      }
      histos.fill(HIST("hnsigmaTPCKaonMinus"), phitrackd1.phid2TPC(), kaonminusd1pt);
      histos.fill(HIST("hPhiMass"), Phid1.M(), Phid1.Pt());
      auto phid1id = phitrackd1.index();
      Phid1.SetXYZM(phitrackd1.phiPx(), phitrackd1.phiPy(), phitrackd1.phiPz(), phitrackd1.phiMass());
      // Phi1kaonplus.SetXYZM(phitrackd1.phid1Px(), phitrackd1.phid1Py(), phitrackd1.phid1Pz(), 0.493);
      // Phi1kaonminus.SetXYZM(phitrackd1.phid2Px(), phitrackd1.phid2Py(), phitrackd1.phid2Pz(), 0.493);
      for (auto phitrackd2 : phitracks) {
        auto phid2id = phitrackd2.index();
        if (phid2id <= phid1id) {
          continue;
        }
        auto kaonplusd2pt = TMath::Sqrt(phitrackd2.phid1Px() * phitrackd2.phid1Px() + phitrackd2.phid1Py() * phitrackd2.phid1Py());
        auto kaonminusd2pt = TMath::Sqrt(phitrackd2.phid2Px() * phitrackd2.phid2Px() + phitrackd2.phid2Py() * phitrackd2.phid2Py());
        if (kaonplusd2pt > maxKaonPt) {
          continue;
        }
        if (kaonminusd2pt > maxKaonPt) {
          continue;
        }
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
        // Phi2kaonplus.SetXYZM(phitrackd2.phid1Px(), phitrackd2.phid1Py(), phitrackd2.phid1Pz(), 0.493);
        // Phi2kaonminus.SetXYZM(phitrackd2.phid2Px(), phitrackd2.phid2Py(), phitrackd2.phid2Pz(), 0.493);

        /*
              // Like
              Phid1like = Phi1kaonplus + Phi2kaonplus;
              Phid2like = Phi1kaonminus + Phi2kaonminus;
              exoticlike = Phid1like + Phid2like;
              auto deltaRlike = TMath::Sqrt(TMath::Power(Phid1like.Phi() - Phid2like.Phi(), 2.0) + TMath::Power(Phid1like.Eta() - Phid2like.Eta(), 2.0));
              auto costhetalike = (Phid1like.Px() * Phid2like.Px() + Phid1like.Py() * Phid2like.Py() + Phid1like.Pz() * Phid2like.Pz()) / (Phid1like.P() * Phid2like.P());
              auto deltamlike = TMath::Sqrt(TMath::Power(Phid1like.M() - 1.0192, 2.0) + TMath::Power(Phid2like.M() - 1.0192, 2.0));
              if (!isDeep) {
                histos.fill(HIST("SEMassLike"), exoticlike.M(), exoticlike.Pt(), deltaRlike, costhetalike, deltamlike, phimult);
              }
              if (isDeep) {
                histos.fill(HIST("SEMassLike"), exoticlike.M(), exoticlike.Pt(), deltaRlike, deepangle(Phid1like, Phid2like), deltamlike, phimult);
              }
        */
        // Unlike
        if (phitrackd2.phiMass() < minPhiMass || phitrackd2.phiMass() > maxPhiMass) {
          continue;
        }
        if (phitrackd1.phiMass() < minPhiMass || phitrackd1.phiMass() > maxPhiMass) {
          continue;
        }
        exotic = Phid1 + Phid2;
        if (exotic.M() < minExoticMass || exotic.M() > maxExoticMass) {
          continue;
        }
        // auto cosThetaStar = getCosTheta(exotic, Phid1);
        // auto kstar = getkstar(Phid1, Phid2);
        auto deltaR = TMath::Sqrt(TMath::Power(Phid1.Phi() - Phid2.Phi(), 2.0) + TMath::Power(Phid1.Eta() - Phid2.Eta(), 2.0));
        // auto costheta = (Phid1.Px() * Phid2.Px() + Phid1.Py() * Phid2.Py() + Phid1.Pz() * Phid2.Pz()) / (Phid1.P() * Phid2.P());
        // auto deltam = TMath::Sqrt(TMath::Power(Phid1.M() - 1.0192, 2.0) + TMath::Power(Phid2.M() - 1.0192, 2.0));
        histos.fill(HIST("hPhiMass2"), Phid1.M(), Phid2.M());
        if (!isDeep) {
          histos.fill(HIST("SEMassUnlike"), exotic.M(), exotic.Pt(), deltaR, Phid1.M(), Phid2.M(), phimult);
        }
        if (isDeep) {
          histos.fill(HIST("SEMassUnlike"), exotic.M(), exotic.Pt(), deltaR, Phid1.M(), Phid2.M(), phimult);
        }
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processSE, "Process Same Event", false);
  void processopti(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    std::vector<ROOT::Math::PtEtaPhiMVector> exoticresonance, phiresonanced1, phiresonanced2;
    std::vector<int> d1trackid = {};
    std::vector<int> d2trackid = {};
    std::vector<int> d3trackid = {};
    std::vector<int> d4trackid = {};
    if (additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }
    int phimult = 0;

    for (auto phitrackd1 : phitracks) {
      if (phitrackd1.phiMass() < minPhiMass || phitrackd1.phiMass() > maxPhiMass) {
        continue;
      }
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
      if (kaonplusd1pt > maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), strategyPID1, kaonminusd1pt)) {
        continue;
      }
      phimult = phimult + 1;
    }
    for (auto phitrackd1 : phitracks) {
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
      if (kaonplusd1pt > maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), strategyPID1, kaonplusd1pt)) {
        continue;
      }
      histos.fill(HIST("hnsigmaTPCTOFKaon"), phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), phitrackd1.phid1TPC(), kaonplusd1pt);
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), strategyPID1, kaonminusd1pt)) {
        continue;
      }
      histos.fill(HIST("hnsigmaTPCKaonMinus"), phitrackd1.phid2TPC(), kaonminusd1pt);
      histos.fill(HIST("hPhiMass"), Phid1.M(), Phid1.Pt());
      auto phid1id = phitrackd1.index();
      Phid1.SetXYZM(phitrackd1.phiPx(), phitrackd1.phiPy(), phitrackd1.phiPz(), phitrackd1.phiMass());
      // Phi1kaonplus.SetXYZM(phitrackd1.phid1Px(), phitrackd1.phid1Py(), phitrackd1.phid1Pz(), 0.493);
      // Phi1kaonminus.SetXYZM(phitrackd1.phid2Px(), phitrackd1.phid2Py(), phitrackd1.phid2Pz(), 0.493);
      for (auto phitrackd2 : phitracks) {
        auto phid2id = phitrackd2.index();
        if (phid2id <= phid1id) {
          continue;
        }
        auto kaonplusd2pt = TMath::Sqrt(phitrackd2.phid1Px() * phitrackd2.phid1Px() + phitrackd2.phid1Py() * phitrackd2.phid1Py());
        auto kaonminusd2pt = TMath::Sqrt(phitrackd2.phid2Px() * phitrackd2.phid2Px() + phitrackd2.phid2Py() * phitrackd2.phid2Py());
        if (kaonplusd2pt > maxKaonPt) {
          continue;
        }
        if (kaonminusd2pt > maxKaonPt) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid1TPC(), phitrackd2.phid1TOF(), phitrackd2.phid1TOFHit(), strategyPID2, kaonplusd2pt)) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid2TPC(), phitrackd2.phid2TOF(), phitrackd2.phid2TOFHit(), strategyPID2, kaonminusd2pt)) {
          continue;
        }
        if ((phitrackd1.phid1Index() == phitrackd2.phid1Index()) || (phitrackd1.phid2Index() == phitrackd2.phid2Index())) {
          continue;
        }
        Phid2.SetXYZM(phitrackd2.phiPx(), phitrackd2.phiPy(), phitrackd2.phiPz(), phitrackd2.phiMass());
        // Phi2kaonplus.SetXYZM(phitrackd2.phid1Px(), phitrackd2.phid1Py(), phitrackd2.phid1Pz(), 0.493);
        // Phi2kaonminus.SetXYZM(phitrackd2.phid2Px(), phitrackd2.phid2Py(), phitrackd2.phid2Pz(), 0.493);
        /*
              // Like
              Phid1like = Phi1kaonplus + Phi2kaonplus;
              Phid2like = Phi1kaonminus + Phi2kaonminus;
              exoticlike = Phid1like + Phid2like;
              auto deltaRlike = TMath::Sqrt(TMath::Power(Phid1like.Phi() - Phid2like.Phi(), 2.0) + TMath::Power(Phid1like.Eta() - Phid2like.Eta(), 2.0));
              auto costhetalike = (Phid1like.Px() * Phid2like.Px() + Phid1like.Py() * Phid2like.Py() + Phid1like.Pz() * Phid2like.Pz()) / (Phid1like.P() * Phid2like.P());
              auto deltamlike = TMath::Sqrt(TMath::Power(Phid1like.M() - 1.0192, 2.0) + TMath::Power(Phid2like.M() - 1.0192, 2.0));
              if (!isDeep) {
                histos.fill(HIST("SEMassLike"), exoticlike.M(), exoticlike.Pt(), deltaRlike, costhetalike, deltamlike, phimult);
              }
              if (isDeep) {
                histos.fill(HIST("SEMassLike"), exoticlike.M(), exoticlike.Pt(), deltaRlike, deepangle(Phid1like, Phid2like), deltamlike, phimult);
              }
        */
        // unlike
        if (phitrackd1.phiMass() < minPhiMass || phitrackd1.phiMass() > maxPhiMass) {
          continue;
        }
        if (phitrackd2.phiMass() < minPhiMass || phitrackd2.phiMass() > maxPhiMass) {
          continue;
        }
        exotic = Phid1 + Phid2;
        if (exotic.M() < minExoticMass || exotic.M() > maxExoticMass) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector temp1(exotic.Pt(), exotic.Eta(), exotic.Phi(), exotic.M());
        ROOT::Math::PtEtaPhiMVector temp2(Phid1.Pt(), Phid1.Eta(), Phid1.Phi(), Phid1.M());
        ROOT::Math::PtEtaPhiMVector temp3(Phid2.Pt(), Phid2.Eta(), Phid2.Phi(), Phid2.M());
        exoticresonance.push_back(temp1);
        phiresonanced1.push_back(temp2);
        phiresonanced2.push_back(temp3);
        d1trackid.push_back(phitrackd1.phid1Index());
        d2trackid.push_back(phitrackd2.phid1Index());
        d3trackid.push_back(phitrackd1.phid2Index());
        d4trackid.push_back(phitrackd2.phid2Index());
      }
    }
    if (exoticresonance.size() == 0) {
      return;
    }
    // LOGF(info, "Total number of exotic: %d", exoticresonance.size());
    if (exoticresonance.size() == 2) {
      for (auto if1 = exoticresonance.begin(); if1 != exoticresonance.end(); ++if1) {
        auto i5 = std::distance(exoticresonance.begin(), if1);
        auto exotic1phi1 = phiresonanced1.at(i5);
        auto exotic1phi2 = phiresonanced2.at(i5);
        auto exotic1 = exoticresonance.at(i5);
        auto deltam1 = TMath::Sqrt(TMath::Power(exotic1phi1.M() - 1.0192, 2.0) + TMath::Power(exotic1phi2.M() - 1.0192, 2.0));
        auto deltaR1 = TMath::Sqrt(TMath::Power(exotic1phi1.Phi() - exotic1phi2.Phi(), 2.0) + TMath::Power(exotic1phi1.Eta() - exotic1phi2.Eta(), 2.0));
        for (auto if2 = if1 + 1; if2 != exoticresonance.end(); ++if2) {
          auto i6 = std::distance(exoticresonance.begin(), if2);
          auto exotic2phi1 = phiresonanced1.at(i6);
          auto exotic2phi2 = phiresonanced2.at(i6);
          auto exotic2 = exoticresonance.at(i6);
          auto deltam2 = TMath::Sqrt(TMath::Power(exotic2phi1.M() - 1.0192, 2.0) + TMath::Power(exotic2phi2.M() - 1.0192, 2.0));
          auto deltaR2 = TMath::Sqrt(TMath::Power(exotic2phi1.Phi() - exotic2phi2.Phi(), 2.0) + TMath::Power(exotic2phi1.Eta() - exotic2phi2.Eta(), 2.0));
          // LOGF(info, "exotic 1 kaon ids %d , %d, %d, %d", d1trackid.at(i5), d2trackid.at(i5), d3trackid.at(i5), d4trackid.at(i5));
          // LOGF(info, "exotic 2 kaon ids %d , %d, %d, %d", d1trackid.at(i6), d2trackid.at(i6), d3trackid.at(i6), d4trackid.at(i6));
          if ((d1trackid.at(i5) == d1trackid.at(i6) || d1trackid.at(i5) == d2trackid.at(i6)) &&
              (d2trackid.at(i5) == d1trackid.at(i6) || d2trackid.at(i5) == d2trackid.at(i6)) &&
              (d3trackid.at(i5) == d3trackid.at(i6) || d3trackid.at(i5) == d4trackid.at(i6)) &&
              (d4trackid.at(i5) == d3trackid.at(i6) || d4trackid.at(i5) == d4trackid.at(i6))) {
            // LOGF(info, "Find Pair %f %f", deltam2, deltam1);
            if (deltam2 < deltam1) {
              histos.fill(HIST("SEMassUnlike"), exotic2.M(), exotic2.Pt(), deltaR2, exotic2phi1.M(), exotic2phi2.M(), phimult);
              // LOGF(info, "Fill exotic Id %d which is pair of Id %d", i6, i5);
            } else {
              histos.fill(HIST("SEMassUnlike"), exotic1.M(), exotic1.Pt(), deltaR1, exotic1phi1.M(), exotic1phi2.M(), phimult);
              // LOGF(info, "Fill exotic Id %d which is pair of Id %d", i6, i5);
            }
          } else {
            histos.fill(HIST("SEMassUnlike"), exotic1.M(), exotic1.Pt(), deltaR1, exotic1phi1.M(), exotic1phi2.M(), phimult);
          }
        }
      }
    } else {
      for (auto if1 = exoticresonance.begin(); if1 != exoticresonance.end(); ++if1) {
        auto i5 = std::distance(exoticresonance.begin(), if1);
        auto exotic1phi1 = phiresonanced1.at(i5);
        auto exotic1phi2 = phiresonanced2.at(i5);
        auto exotic1 = exoticresonance.at(i5);
        // auto deltam1 = TMath::Sqrt(TMath::Power(exotic1phi1.M() - 1.0192, 2.0) + TMath::Power(exotic1phi2.M() - 1.0192, 2.0));
        auto deltaR1 = TMath::Sqrt(TMath::Power(exotic1phi1.Phi() - exotic1phi2.Phi(), 2.0) + TMath::Power(exotic1phi1.Eta() - exotic1phi2.Eta(), 2.0));
        histos.fill(HIST("SEMassUnlike"), exotic1.M(), exotic1.Pt(), deltaR1, exotic1phi1.M(), exotic1phi2.M(), phimult);
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processopti, "Process Optimized same event", false);

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
        if (kaonplusd1pt > maxKaonPt) {
          continue;
        }
        if (kaonminusd1pt > maxKaonPt) {
          continue;
        }
        if (kaonplusd2pt > maxKaonPt) {
          continue;
        }
        if (kaonminusd2pt > maxKaonPt) {
          continue;
        }
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
        auto deltaR = TMath::Sqrt(TMath::Power(Phid1.Phi() - Phid2.Phi(), 2.0) + TMath::Power(Phid1.Eta() - Phid2.Eta(), 2.0));
        // auto costheta = (Phid1.Px() * Phid2.Px() + Phid1.Py() * Phid2.Py() + Phid1.Pz() * Phid2.Pz()) / (Phid1.P() * Phid2.P());
        // auto deltam = TMath::Sqrt(TMath::Power(Phid1.M() - 1.0192, 2.0) + TMath::Power(Phid2.M() - 1.0192, 2.0));
        if (!isDeep) {
          histos.fill(HIST("MEMassUnlike"), exotic.M(), exotic.Pt(), deltaR, Phid1.M(), Phid2.M());
        }
        if (isDeep) {
          histos.fill(HIST("MEMassUnlike"), exotic.M(), exotic.Pt(), deltaR, Phid1.M(), Phid2.M());
        }
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processMixedEvent, "Process EventMixing for combinatorial background", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<doublephimeson>(cfgc)}; }
